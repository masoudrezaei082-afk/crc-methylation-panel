"""CpG-gene-TF multiplex network (STRING, JASPAR) and aggregated multiplex PageRank."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

LAYERS = ("cpg", "ppi", "tf")


def _request(method: str, url: str, attempts: int = 10, **kw):
    """HTTP call that survives DNS failures, timeouts and 5xx (backoff up to ~20 min total)."""
    import requests

    for i in range(attempts):
        try:
            r = requests.request(method, url, **kw)
            if r.status_code < 500 and r.status_code != 429:
                r.raise_for_status()
                return r
            err = f"HTTP {r.status_code}"
        except requests.exceptions.HTTPError:
            raise
        except requests.exceptions.RequestException as exc:
            err = repr(exc)
        if i == attempts - 1:
            raise RuntimeError(f"{method} {url} failed after {attempts} attempts: {err}")
        wait = min(300, 5 * 2 ** i)
        print(f"  {method} {url[:80]}: {err[:120]}; retry in {wait}s", flush=True)
        time.sleep(wait)


# --------------------------------------------------------------------------- STRING
def string_edges(genes, cfg: dict, cache_dir: Path, required_score: int | None = None) -> pd.DataFrame:
    """Undirected PPI edges among `genes` (STRING combined score / 1000)."""
    genes = sorted(set(genes))
    score = required_score or cfg["network"]["string_required_score"]
    key = hashlib.md5(("|".join(genes) + str(score)).encode()).hexdigest()[:16]
    cache = Path(cache_dir) / f"string_{key}.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    api = cfg["urls"]["string_api"]
    ids = _request("POST", f"{api}/tsv/get_string_ids",
                   data={"identifiers": "\r".join(genes), "species": 9606, "limit": 1,
                         "echo_query": 1, "caller_identity": "crcmeth_repro"},
                   timeout=300)
    idmap = pd.read_csv(pd.io.common.StringIO(ids.text), sep="\t")
    to_query = dict(zip(idmap["preferredName"], idmap["queryItem"]))
    time.sleep(1)
    r = _request("POST", f"{api}/tsv/network",
                 data={"identifiers": "\r".join(idmap["stringId"]), "species": 9606,
                       "required_score": score, "caller_identity": "crcmeth_repro"},
                 timeout=300)
    net = pd.read_csv(pd.io.common.StringIO(r.text), sep="\t")
    out = pd.DataFrame({"a": net["preferredName_A"].map(to_query),
                        "b": net["preferredName_B"].map(to_query),
                        "weight": net["score"].astype(float)}).dropna()
    out = out[out["a"] != out["b"]].drop_duplicates()
    out.to_csv(cache, sep="\t", index=False)
    return out


# --------------------------------------------------------------------------- JASPAR
def promoter_sequence(chrom: str, tss: int, strand: str, cfg: dict, cache_dir: Path) -> str:
    """hg38 promoter (-up/+down around TSS) in transcript orientation, UCSC REST."""
    up, down = cfg["network"]["promoter_upstream"], cfg["network"]["promoter_downstream"]
    start, end = (tss - up, tss + down) if strand == "+" else (tss - down, tss + up)
    cache = Path(cache_dir) / "promoters" / f"{chrom}_{start}_{end}.txt"
    if cache.exists():
        seq = cache.read_text()
    else:
        url = f"{cfg['urls']['ucsc_api']}?genome=hg38;chrom={chrom};start={start - 1};end={end}"
        r = _request("GET", url, timeout=60)
        seq = r.json()["dna"].upper()
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(seq)
        time.sleep(0.3)                     # UCSC asks for <= ~1 request/second bursts
    if strand == "-":
        from Bio.Seq import Seq
        seq = str(Seq(seq).reverse_complement())
    return seq


def load_pssms(cfg: dict, cache_dir: Path):
    """JASPAR CORE vertebrate PSSMs with a p-value score threshold each."""
    import pickle

    from pyjaspar import jaspardb

    n = cfg["network"]
    cache = Path(cache_dir) / f"pssm_{n['jaspar_release']}_{n['motif_pvalue']:g}.pkl"
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    jdb = jaspardb(release=n["jaspar_release"])
    motifs = jdb.fetch_motifs(collection=n["jaspar_collection"],
                              tax_group=[n["jaspar_tax_group"]])
    out = []
    for m in motifs:
        pssm = m.counts.normalize(pseudocounts=0.5).log_odds()
        thr = pssm.distribution(precision=10 ** 3).threshold_fpr(n["motif_pvalue"])
        tfs = sorted({t.upper() for t in m.name.split("::")})
        out.append((m.matrix_id, tfs, pssm, float(thr)))
    cache.write_bytes(pickle.dumps(out))
    return out


def tf_edges(genes, probemap: pd.DataFrame, cfg: dict, cache_dir: Path) -> pd.DataFrame:
    """TF -> gene edges: >=1 motif hit (p < motif_pvalue) in the promoter, weight 1."""
    from Bio.Seq import Seq

    cache = Path(cache_dir) / "tf_hits"
    cache.mkdir(parents=True, exist_ok=True)
    pssms = None
    pm = probemap.drop_duplicates("gene").set_index("gene")
    rows = []
    for g in sorted(set(genes)):
        f = cache / f"{g}.json"
        if f.exists():
            hits = json.loads(f.read_text())
        else:
            if g not in pm.index:
                continue
            pssms = pssms or load_pssms(cfg, cache_dir)
            rec = pm.loc[g]
            seq = Seq(promoter_sequence(rec["chrom"], int(rec["tss"]), rec["strand"], cfg, cache_dir))
            hits = sorted({tf for _, tfs, pssm, thr in pssms
                           if np.nanmax(np.concatenate([np.atleast_1d(pssm.calculate(seq)),
                                                        np.atleast_1d(pssm.reverse_complement().calculate(seq))])) >= thr
                           for tf in tfs})
            f.write_text(json.dumps(hits))
        rows.extend((tf, g) for tf in hits)
    return pd.DataFrame(rows, columns=["tf", "gene"]).assign(weight=1.0)


# --------------------------------------------------------------------------- graph
def build_graph(cpg_gene: pd.DataFrame, ppi: pd.DataFrame, tf: pd.DataFrame,
                weights: dict) -> nx.DiGraph:
    """Aggregated multiplex A* = sum_l beta_l A^(l) (paper Eq. 2).

    Node ids carry a type prefix so a TF and the gene of the same symbol stay
    distinct nodes, as in the manuscript's counts (genes, CpGs, TFs separately)."""
    G = nx.DiGraph()
    genes = set(cpg_gene["gene"])
    G.add_nodes_from((f"G:{g}" for g in genes), kind="gene")

    def add(u, v, w):
        if G.has_edge(u, v):
            G[u][v]["weight"] += w
        else:
            G.add_edge(u, v, weight=w)

    for r in cpg_gene.itertuples():
        G.add_node(f"C:{r.probe}", kind="cpg")
        add(f"C:{r.probe}", f"G:{r.gene}", weights["cpg"] * r.weight)
    for r in ppi.itertuples():
        if r.a in genes and r.b in genes:
            add(f"G:{r.a}", f"G:{r.b}", weights["ppi"] * r.weight)
            add(f"G:{r.b}", f"G:{r.a}", weights["ppi"] * r.weight)
    for r in tf.itertuples():
        if r.gene in genes:
            G.add_node(f"T:{r.tf}", kind="tf")
            add(f"T:{r.tf}", f"G:{r.gene}", weights["tf"] * r.weight)
    return G


def pagerank(G: nx.DiGraph, alpha: float = 0.85, tol: float = 1e-6) -> pd.Series:
    """Power iteration; dangling mass to the uniform vector (networkx convention).

    Converges when the L1 change is below tol (the paper's 1e-6 criterion)."""
    nodes = list(G.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    N = len(nodes)
    W = np.zeros((N, N))
    for u, v, d in G.edges(data=True):
        W[idx[v], idx[u]] += d["weight"]
    out = W.sum(axis=0)
    dangling = out == 0
    W[:, ~dangling] /= out[~dangling]
    v = np.full(N, 1.0 / N)
    p = v.copy()
    for _ in range(10000):
        new = alpha * (W @ p + p[dangling].sum() / N) + (1 - alpha) * v
        if np.abs(new - p).sum() < tol:
            p = new
            break
        p = new
    return pd.Series(p / p.sum(), index=nodes)


def gene_ranking(G: nx.DiGraph, cfg: dict) -> pd.DataFrame:
    n = cfg["network"]
    pr = pagerank(G, n["alpha"], n["tol"])
    nx_pr = pd.Series(nx.pagerank(G, alpha=n["alpha"], tol=n["tol"] / len(G), max_iter=1000, weight="weight"))
    genes = [x for x in pr.index if x.startswith("G:")]
    df = pd.DataFrame({"pagerank": pr[genes], "pagerank_networkx": nx_pr[genes]})
    df.index = [g[2:] for g in genes]
    df["in_cpg"] = [sum(1 for u in G.predecessors(f"G:{g}") if u.startswith("C:")) for g in df.index]
    df["in_tf"] = [sum(1 for u in G.predecessors(f"G:{g}") if u.startswith("T:")) for g in df.index]
    df["ppi_degree"] = [sum(1 for u in G.successors(f"G:{g}") if u.startswith("G:")) for g in df.index]
    df = df.sort_values("pagerank", ascending=False)
    df["rank"] = np.arange(1, len(df) + 1)
    df["tie_group_size"] = df.groupby(df["pagerank"].round(10))["pagerank"].transform("size")
    return df


def beta_sweep(cpg_gene, ppi, tf, cfg: dict, reference_top) -> pd.DataFrame:
    """Top-k overlap and Kendall tau vs reference over the layer-weight simplex."""
    from scipy.stats import kendalltau

    step = cfg["network"]["beta_grid_step"]
    k = cfg["network"]["top_k_genes"]
    grid = np.round(np.arange(0, 1 + 1e-9, step), 6)
    rows = []
    base = None
    for b1 in grid:
        for b2 in grid:
            b3 = round(1 - b1 - b2, 6)
            if b3 < -1e-9 or b1 == 0:          # CpG layer must exist to reach genes from CpGs
                continue
            w = {"cpg": b1, "ppi": b2, "tf": max(b3, 0.0)}
            rk = gene_ranking(build_graph(cpg_gene, ppi, tf, w), cfg)
            if base is None:
                base = rk.index
            top = list(rk.index[:k])
            tau = kendalltau(rk["pagerank"].reindex(base).values,
                             np.arange(len(base))[::-1]).statistic
            rows.append({**w, "top_k": ";".join(top),
                         "overlap_with_reference": len(set(top) & set(reference_top)),
                         "max_tie_in_top_k": int(rk["tie_group_size"].iloc[:k].max()),
                         "kendall_vs_first": tau})
    return pd.DataFrame(rows)


def edge_bootstrap(cpg_gene, ppi, tf, cfg: dict, rng: np.random.Generator) -> pd.DataFrame:
    """Frequency with which each gene stays in the top-k after random edge removal."""
    n = cfg["network"]
    k, keep = n["top_k_genes"], n["edge_keep_frac"]
    counts: dict[str, int] = {}
    for _ in range(n["n_bootstrap"]):
        sub = [df.loc[rng.random(len(df)) < keep] for df in (cpg_gene, ppi, tf)]
        rk = gene_ranking(build_graph(*sub, n["layer_weights"]), cfg)
        for g in rk.index[:k]:
            counts[g] = counts.get(g, 0) + 1
    return (pd.Series(counts, name="top_k_frequency") / n["n_bootstrap"]).sort_values(ascending=False).to_frame()
