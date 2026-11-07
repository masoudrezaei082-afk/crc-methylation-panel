#!/usr/bin/env python
"""End-to-end, cached pipeline for the CRC methylation-marker analysis.

    python run_all.py                      # every stage, cached
    python run_all.py --stage classify     # one stage (its inputs must be cached)
    python run_all.py --from dmc --force   # recompute dmc and everything after it

Stage outputs go to <workdir>/cache; human-readable results to <resultsdir>.
"""

from __future__ import annotations

import argparse
import hashlib
import pickle
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from crcmeth import annot as A          # noqa: E402
from crcmeth import blood as B          # noqa: E402
from crcmeth import classify as C       # noqa: E402
from crcmeth import data as D           # noqa: E402
from crcmeth import deg as G            # noqa: E402
from crcmeth import dmc as M            # noqa: E402
from crcmeth import integrate as I      # noqa: E402
from crcmeth import network as N        # noqa: E402
from crcmeth import report as R         # noqa: E402
from crcmeth import specificity as S    # noqa: E402

STAGES = ["annot", "coad", "pancan", "dmc", "rna", "integrate", "specificity",
          "network", "geo", "classify", "blood", "report"]


# --------------------------------------------------------------------------- cache helpers
class Ctx:
    def __init__(self, cfg: dict, force: set[str]):
        self.cfg = cfg
        self.cache = Path(cfg["workdir"]) / "cache"
        self.cache.mkdir(parents=True, exist_ok=True)
        self.res = Path(cfg["resultsdir"])
        for sub in ("tables", "figures", "logs"):
            (self.res / sub).mkdir(parents=True, exist_ok=True)
        self.force = force

    def path(self, name: str) -> Path:
        return self.cache / name

    def done(self, stage: str) -> bool:
        return stage not in self.force and self.path(f"{stage}.done").exists()

    def mark(self, stage: str, info: dict) -> None:
        D.write_json(info, self.path(f"{stage}.done"))
        D.write_json(info, self.res / "logs" / f"stage_{stage}.json")

    def dump(self, name: str, obj) -> None:
        self.path(name).write_bytes(pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL))

    def load(self, name: str):
        return pickle.loads(self.path(name).read_bytes())

    def cached(self, stage: str, name: str, fn):
        """Checkpoint inside a long stage: reuse `name` if a crashed run already wrote it."""
        if stage not in self.force and self.path(name).exists():
            return self.load(name)
        obj = fn()
        self.dump(name, obj)
        return obj

    def load_info(self, stage: str) -> dict:
        import json
        return json.loads(self.path(f"{stage}.done").read_text())

    def table(self, df: pd.DataFrame, name: str, index: bool = True) -> None:
        df.to_csv(self.res / "tables" / name, index=index)


def log(msg: str) -> None:
    print(time.strftime("[%H:%M:%S] ") + msg, flush=True)


def keep_raw(ctx: Ctx, path: Path) -> None:
    if not ctx.cfg["keep_raw"]:
        path.unlink(missing_ok=True)


# --------------------------------------------------------------------------- stages
def stage_annot(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    ann = A.load_gpl13534(D.download(cfg, cfg["urls"]["gpl13534"]))
    mask = A.load_zhou_mask(D.download(cfg, cfg["urls"]["zhou_mask"]))
    pm = A.gene_probemap(D.download(cfg, cfg["urls"]["gene_probemap"], "gencode_v36_probemap.tsv"))
    ctx.dump("annot.pkl", ann)
    ctx.dump("mask.pkl", mask)
    ctx.dump("probemap.pkl", pm)
    sept9 = A.genes_in(ann, cfg["promoter_groups"]["primary"])
    sept9 = sorted(sept9.index[sept9.map(lambda g: "SEPT9" in g)])
    ctx.dump("sept9.pkl", sept9)
    return {"n_probes": len(ann), "n_masked": len(mask), "n_genes_probemap": len(pm),
            "n_sept9_promoter_probes": len(sept9)}


def stage_coad(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    path = D.download(cfg, D.gdc_url(cfg, "COAD", "methylation450"))
    cols = D.header_columns(path)
    tumors = D.select_samples(cols, D.TUMOR_CODE)
    normals = D.select_samples(cols, D.NORMAL_CODE)
    log(f"COAD 450K: {len(cols) - 1} columns -> {len(tumors)} tumours, {len(normals)} normals")
    beta = D.read_matrix(path, tumors + normals, cfg["dmc"]["chunk_rows"])
    D.save_frame(beta, ctx.path("coad_beta_raw.parquet"))
    keep_raw(ctx, path)
    return {"n_columns": len(cols) - 1, "n_tumor": len(tumors), "n_normal": len(normals),
            "tumors": tumors, "normals": normals, "n_probes_raw": int(beta.shape[0])}


def stage_pancan(ctx: Ctx) -> dict:
    """Stream each other cancer once: genome-wide stats, kept rows, Design A normals."""
    cfg = ctx.cfg
    coad = ctx.load_info("coad")
    beta = D.load_frame(ctx.path("coad_beta_raw.parquet"))
    delta = beta[coad["tumors"]].mean(axis=1) - beta[coad["normals"]].mean(axis=1)
    keep = set(delta.index[delta > 0.15]) | set(ctx.load("sept9.pkl"))
    del beta
    cancers = list(dict.fromkeys(cfg["screen_cancers"] + cfg["design_a_normal_sources"]))
    meta = {}
    out_dir = ctx.path("pancan")
    out_dir.mkdir(exist_ok=True)
    for i, cancer in enumerate(cancers):
        if (out_dir / f"{cancer}.meta.pkl").exists() and "pancan" not in ctx.force:
            meta[cancer] = pickle.loads((out_dir / f"{cancer}.meta.pkl").read_bytes())
            continue
        t0 = time.time()
        path = D.download(cfg, D.gdc_url(cfg, cancer, "methylation450"))
        n_ext = cfg["design_a_per_source"] if cancer in cfg["design_a_normal_sources"] else 0
        st, kept, ext, m = S.pancancer_stats(path, cfg, keep, n_ext, D.rng(cfg, 100 + i))
        D.save_frame(st, out_dir / f"{cancer}.stats.parquet")
        D.save_frame(kept, out_dir / f"{cancer}.kept.parquet")
        if ext is not None:
            D.save_frame(ext, out_dir / f"{cancer}.normals.parquet")
        (out_dir / f"{cancer}.meta.pkl").write_bytes(pickle.dumps(m))
        meta[cancer] = m
        keep_raw(ctx, path)
        log(f"  {cancer}: {m['n_tumor']} T / {m['n_normal']} N, {time.time() - t0:.0f}s")
    rows = [{"cancer": c, "n_tumor": m["n_tumor"], "n_normal": m["n_normal"],
             "in_screen": c in cfg["screen_cancers"],
             "design_a_normals": ";".join(m["extracted_normals"])} for c, m in meta.items()]
    ctx.table(pd.DataFrame(rows), "pancancer_cohorts.csv", index=False)
    return {"cancers": meta, "n_keep_probes": len(keep)}


def _qc_impute(ctx: Ctx, beta: pd.DataFrame) -> pd.DataFrame:
    """Discovery matrix: all QC and imputation fitted on every sample given."""
    b = M.probe_qc(beta, ctx.load("annot.pkl"), ctx.load("mask.pkl"), ctx.cfg, missing=False)
    return M.filter_impute(b, ctx.cfg)


def _model_matrix(ctx: Ctx, samples) -> pd.DataFrame:
    """Evaluation matrix: COAD betas after the sample-independent probe filters only. The
    missingness filter and imputation are left to each training fold (nested funnel and
    model pipeline), so no test sample shapes the values a model is trained on."""
    raw = D.load_frame(ctx.path("coad_beta_raw.parquet"))[samples]
    return M.probe_qc(raw, ctx.load("annot.pkl"), ctx.load("mask.pkl"), ctx.cfg, missing=False)


def _patients(samples) -> np.ndarray:
    return np.array([D.patient_of(s) for s in samples])


def stage_dmc(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    coad = ctx.load_info("coad")
    pan = ctx.load_info("pancan")["cancers"]
    raw = D.load_frame(ctx.path("coad_beta_raw.parquet"))
    extra = [D.load_frame(ctx.path("pancan") / f"{c}.normals.parquet")
             for c in cfg["design_a_normal_sources"]]
    extra = pd.concat(extra, axis=1).reindex(raw.index)
    ext_cols = list(extra.columns)
    assert len(ext_cols) == len(cfg["design_a_normal_sources"]) * cfg["design_a_per_source"]
    designs = {"B": coad["normals"], "A": coad["normals"] + ext_cols}
    full = pd.concat([raw, extra], axis=1)
    del raw, extra
    info, tabs = {}, {}
    for name, normals in designs.items():
        b = _qc_impute(ctx, full[coad["tumors"] + normals])
        dt = M.dmc_table(b, coad["tumors"], normals, cfg)
        D.save_frame(dt, ctx.path(f"dmc_{name}.parquet"))
        if name == cfg["designs"]["primary"]:
            D.save_frame(b, ctx.path("coad_beta.parquet"))      # model matrix
        tabs[name] = dt
        vc = dt["status"].value_counts()
        info[name] = {"n_tumor": len(coad["tumors"]), "n_normal": len(normals),
                      "n_probes_tested": int(len(dt)), "n_hyper": int(vc.get("hyper", 0)),
                      "n_hypo": int(vc.get("hypo", 0))}
        info[name]["n_dmc"] = info[name]["n_hyper"] + info[name]["n_hypo"]
        info[name]["regions"] = I.region_summary(dt, ctx.load("annot.pkl"),
                                                 cfg["promoter_groups"]["primary"])
        log(f"  design {name}: {info[name]['n_dmc']} DMCs "
            f"({info[name]['n_hyper']} hyper / {info[name]['n_hypo']} hypo)")
    common = tabs["A"].index.intersection(tabs["B"].index)
    for s in ("hyper", "hypo"):
        a = set(common[tabs["A"].loc[common, "status"] == s])
        b = set(common[tabs["B"].loc[common, "status"] == s])
        info[f"jaccard_{s}"] = len(a & b) / max(1, len(a | b))
        info[f"only_A_{s}"], info[f"only_B_{s}"] = len(a - b), len(b - a)
    info["extracted_normals"] = {c: pan[c]["extracted_normals"] for c in cfg["design_a_normal_sources"]}
    return info


def stage_rna(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    path = D.download(cfg, D.gdc_url(cfg, "COAD", "star_counts"))
    cols = D.header_columns(path)
    tumors = D.select_samples(cols, D.TUMOR_CODE)
    normals = D.select_samples(cols, D.NORMAL_CODE)
    log2c = pd.read_csv(path, sep="\t", index_col=0, usecols=[cols[0]] + tumors + normals)
    counts = G.counts_from_xena(log2c)
    counts = counts.loc[~counts.index.str.startswith("__")]
    D.save_frame(counts, ctx.path("rna_counts.parquet"))
    cond = pd.Series(["tumor"] * len(tumors) + ["normal"] * len(normals), index=tumors + normals)
    filt = G.filter_genes(counts, cfg)
    log(f"  DESeq2 on {filt.shape[0]} genes, {len(tumors)} T / {len(normals)} N")
    res = G.run_deseq2(filt, cond, cfg)
    res = G.to_symbols(res, ctx.load("probemap.pkl"))
    res["status"] = G.label_deg(res, cfg["deg"]["padj"], cfg["deg"]["min_abs_log2fc"])
    D.save_frame(res, ctx.path("deg.parquet"))
    ctx.table(res.loc[res["status"] != "ns"], "deg_significant.csv")
    keep_raw(ctx, path)
    vc = res["status"].value_counts()
    sens = {}
    for lfc in (1.0, 1.5, 2.0):
        lab = G.label_deg(res, cfg["deg"]["padj"], lfc).value_counts()
        sens[str(lfc)] = {"up": int(lab.get("up", 0)), "down": int(lab.get("down", 0))}
    return {"n_tumor": len(tumors), "n_normal": len(normals), "tumors": tumors,
            "normals": normals, "n_genes_tested": int(len(res)),
            "n_up": int(vc.get("up", 0)), "n_down": int(vc.get("down", 0)),
            "lfc_sensitivity": sens}


def _integrate(ctx: Ctx, dmc: pd.DataFrame, deg_status: pd.Series, groups) -> dict:
    ann = ctx.load("annot.pkl")
    dmg = I.promoter_dmg(dmc, ann, groups)
    q = I.quadrants(dmg, deg_status)
    cands = I.candidate_cpgs(dmc, ann, q["hyper_down"], groups)
    return {"dmg": dmg, "quadrants": q, "candidates": cands}


def stage_integrate(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    deg = D.load_frame(ctx.path("deg.parquet"))
    info = {}
    for design in ("A", "B"):
        dmc = D.load_frame(ctx.path(f"dmc_{design}.parquet"))
        for pdef in ("primary", "extended"):
            r = _integrate(ctx, dmc, deg["status"], cfg["promoter_groups"][pdef])
            ctx.dump(f"integrate_{design}_{pdef}.pkl", r)
            q = r["quadrants"]
            info[f"{design}_{pdef}"] = {
                "hyper_genes": q["hyper_genes"], "hypo_genes": q["hypo_genes"],
                **{k: len(q[k]) for k in ("hyper_up", "hyper_down", "hypo_up", "hypo_down")},
                "candidate_cpgs": int(len(r["candidates"]))}
            ctx.table(r["candidates"], f"candidates_{design}_{pdef}.csv")
    prim = ctx.load(f"integrate_{cfg['designs']['primary']}_primary.pkl")
    # methylation-expression correlation over patients profiled by both assays
    beta = D.load_frame(ctx.path("coad_beta.parquet"))
    counts = D.load_frame(ctx.path("rna_counts.parquet"))
    pm = ctx.load("probemap.pkl")
    logcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
    logcpm = logcpm.join(pm["gene"]).groupby("gene").max()
    key = lambda c: "-".join(c.split("-")[:4])[:15]        # patient + sample type
    b = beta.loc[beta.index.intersection(prim["candidates"].index)]
    b.columns = [key(c) for c in b.columns]
    logcpm.columns = [key(c) for c in logcpm.columns]
    cg = prim["candidates"]["genes"].str.split(";").explode()
    rho = I.matched_spearman(b, logcpm, cg)
    ctx.table(rho, "meth_expr_spearman.csv", index=False)
    info["spearman_n_pairs"] = int(rho["n"].iloc[0]) if len(rho) else 0
    info["spearman_median_rho"] = float(rho["rho"].median()) if len(rho) else None
    info["spearman_frac_negative_sig"] = float(((rho["rho"] < 0) & (rho["p"] < 0.05)).mean()) if len(rho) else None
    try:
        enr = I.enrichment(prim["quadrants"]["hyper_down"], deg.index,
                           ["GO_Biological_Process_2023", "KEGG_2021_Human"])
        ctx.table(enr, "enrichment_hyper_down.csv", index=False)
        info["enrichment_terms_fdr05"] = int((enr["Adjusted P-value"] < 0.05).sum())
    except Exception as exc:                                  # network-dependent, non-fatal
        info["enrichment_error"] = repr(exc)
    return info


def _screen(ctx: Ctx, cands: pd.DataFrame, cfg: dict | None = None) -> pd.DataFrame:
    cfg = cfg or ctx.cfg
    pan = {c: D.load_frame(ctx.path("pancan") / f"{c}.stats.parquet") for c in cfg["screen_cancers"]}
    meta = {c: m for c, m in ctx.load_info("pancan")["cancers"].items() if c in cfg["screen_cancers"]}
    return S.screen(cands.index, pan, meta, cfg)


def stage_specificity(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    info = {}
    for design in ("A", "B"):
        for pdef in ("primary", "extended"):
            r = ctx.load(f"integrate_{design}_{pdef}.pkl")
            scr = _screen(ctx, r["candidates"])
            alt = _screen(ctx, r["candidates"],
                          {**cfg, "specificity": {**cfg["specificity"], "rule": "dmc_or_high"}})
            scr["specific_dmc_or_high"] = alt["specific"]
            ctx.dump(f"screen_{design}_{pdef}.pkl", scr)
            ctx.table(scr, f"specificity_{design}_{pdef}.csv")
            info[f"{design}_{pdef}"] = {
                "candidates": int(len(scr)), "specific": int(scr["specific"].sum()),
                "specific_dmc_or_high": int(scr["specific_dmc_or_high"].sum()),
                "excluded_by": {c.replace("_excluded", ""): int(scr[c].sum())
                                for c in scr.columns if c.endswith("_excluded")},
                "underpowered": scr.attrs.get("underpowered", [])}
    return info


def _specific(ctx: Ctx, design: str) -> pd.DataFrame:
    """Candidate CpGs (promoter definition of the panels) that pass the specificity screen."""
    pdef = ctx.cfg["classify"]["panel_promoter_groups"]
    r = ctx.load(f"integrate_{design}_{pdef}.pkl")
    scr = ctx.load(f"screen_{design}_{pdef}.pkl")
    specific = r["candidates"].loc[scr.index[scr["specific"]]]
    assert len(specific), f"design {design}: no candidate CpG survives the specificity screen"
    return specific


def _comparator(ctx: Ctx) -> list[str]:
    """The network-prioritized comparator panel, derived by stage_network."""
    return list(ctx.load_info("network")["panel"])


def _network_inputs(ctx: Ctx, specific: pd.DataFrame, required_score=None):
    cg = specific.assign(gene=specific["genes"].str.split(";")).explode("gene")
    cpg_gene = pd.DataFrame({"probe": cg.index, "gene": cg["gene"].values,
                             "weight": cg["delta"].abs().values})
    genes = sorted(cpg_gene["gene"].unique())
    netc = ctx.path("net_cache")
    netc.mkdir(exist_ok=True)
    ppi = N.string_edges(genes, ctx.cfg, netc, required_score)
    tf = N.tf_edges(genes, ctx.load("probemap.pkl"), ctx.cfg, netc)
    return cpg_gene, ppi, tf


def stage_network(ctx: Ctx) -> dict:
    """CpG-gene-TF multiplex network of the primary design's colon-specific CpGs, ranked by
    aggregated multiplex PageRank; the comparator panel follows from the ranking (C.network_panel)."""
    cfg = ctx.cfg
    n = cfg["network"]
    specific = _specific(ctx, cfg["designs"]["primary"])
    cpg_gene, ppi, tf = _network_inputs(ctx, specific)
    G_ = N.build_graph(cpg_gene, ppi, tf, n["layer_weights"])
    rk = N.gene_ranking(G_, cfg)
    panel = C.network_panel(rk, specific, cfg["classify"]["panel_size"])
    top = list(rk.index[:n["top_k_genes"]])
    ann = ctx.load("annot.pkl").loc[panel].copy()
    ann["gene_group"] = ann["gene_group"].map(lambda gg: ";".join(f"{g}:{r}" for g, r in gg))
    ann["zhou_masked"] = ann.index.isin(ctx.load("mask.pkl"))
    ann["pagerank_rank"] = [";".join(f"{g}:{int(rk.loc[g, 'rank'])}" for g in specific.at[p, "genes"].split(";")
                                     if g in rk.index) for p in ann.index]
    ctx.table(ann, "network_prioritized_panel_annotation.csv")
    ctx.table(rk, "pagerank_network.csv")
    log("  sensitivity: layer-weight sweep, STRING 700, edge bootstrap")
    sweep = N.beta_sweep(cpg_gene, ppi, tf, cfg, top)
    ctx.table(sweep, "pagerank_beta_sweep.csv", index=False)
    ppi7 = N.string_edges(sorted(cpg_gene["gene"].unique()), cfg, ctx.path("net_cache"), 700)
    rk7 = N.gene_ranking(N.build_graph(cpg_gene, ppi7, tf, n["layer_weights"]), cfg)
    ctx.table(rk7, "pagerank_string700.csv")
    boot = N.edge_bootstrap(cpg_gene, ppi, tf, cfg, D.rng(cfg, 7))
    ctx.table(boot, "pagerank_edge_bootstrap.csv")
    kinds = pd.Series(dict(G_.nodes(data="kind"))).value_counts()
    try:                                     # the release the live API served
        string_release = N._request("GET", f"{cfg['urls']['string_api']}/json/version", timeout=60).json()
    except Exception as exc:                 # recorded, not fatal: the edge cache is what was used
        string_release = repr(exc)
    return {"panel": panel, "top10": top, "top10_pagerank": rk["pagerank"].iloc[:len(top)].round(8).tolist(),
            "max_tie_in_top10": int(rk["tie_group_size"].iloc[:len(top)].max()),
            "n_specific_cpgs": int(len(specific)), "n_genes": int(kinds.get("gene", 0)),
            "n_cpgs": int(kinds.get("cpg", 0)), "n_tfs": int(kinds.get("tf", 0)),
            "n_edges": G_.number_of_edges(), "n_ppi_edges": int(len(ppi)), "n_tf_edges": int(len(tf)),
            "max_abs_diff_vs_networkx": float((rk["pagerank"] - rk["pagerank_networkx"]).abs().max()),
            "sweep_min_top_overlap": int(sweep["overlap_with_reference"].min()),
            "sweep_median_top_overlap": float(sweep["overlap_with_reference"].median()),
            "string700_top10": list(rk7.index[:len(top)]),
            "string700_top10_overlap": len(set(rk7.index[:len(top)]) & set(top)),
            "bootstrap_top10_freq": boot["top_k_frequency"].reindex(top).fillna(0).round(3).tolist(),
            "panel_masked": int(ann["zhou_masked"].sum()),
            "parameters": {**n, "string_release_served": string_release}}


def _read_gse131013(path: Path, probes) -> pd.DataFrame:
    head = D.header_columns(path)
    vals = [c for c in head[1:] if not c.endswith("Pval")]
    out = []
    for chunk in pd.read_csv(path, sep="\t", index_col=0, usecols=[head[0]] + vals, chunksize=50000):
        out.append(chunk.loc[chunk.index.isin(probes)])
    return pd.concat(out).astype(np.float32)


def stage_geo(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    prim = cfg["designs"]["primary"]
    probes = (set(_comparator(ctx)) | set(ctx.load("sept9.pkl"))
              | set(_specific(ctx, "A").index) | set(_specific(ctx, "B").index))
    info = {}
    # GSE48684: betas inline in the series matrix
    md, tab = D.parse_series_matrix(D.download(cfg, cfg["urls"]["gse48684_matrix"]))
    status = md["disease status"].str.lower()
    lab = pd.Series(np.select([status.str.startswith("cancer"), status.str.startswith("normal")],
                              [1.0, 0.0], np.nan), index=md.index)
    grp = status.rename("group")
    b48 = tab.loc[tab.index.isin(probes)]
    D.save_frame(b48, ctx.path("gse48684_beta.parquet"))
    ctx.dump("gse48684_labels.pkl", pd.DataFrame({"y": lab, "group": grp}))
    info["GSE48684"] = {"n_samples": int(len(md)), "groups": grp.value_counts().to_dict(),
                        "beta_range": [float(np.nanmin(tab.values)), float(np.nanmax(tab.values))],
                        "missing_network_prioritized": sorted(set(_comparator(ctx)) - set(b48.index))}
    del tab
    # GSE131013: betas only in the supplementary normalized matrix
    md2, _ = D.parse_series_matrix(D.download(cfg, cfg["urls"]["gse131013_matrix"]), with_table=False)
    b131 = _read_gse131013(D.download(cfg, cfg["urls"]["gse131013_betas"]), probes)
    # <patient>_<T|N|M>, one sample (Z2015) carries an extra replicate suffix "_0"
    suffix = pd.Series(b131.columns.str.extract(r"_([TNM])(?:_\d+)?$")[0].values, index=b131.columns)
    assert suffix.notna().all(), list(b131.columns[suffix.isna().values])
    grp2 = suffix.map({"T": "tumor", "N": "adjacent_normal", "M": "healthy_mucosa"})
    lab2 = grp2.map({"tumor": 1.0, "adjacent_normal": 0.0, "healthy_mucosa": 0.0})
    D.save_frame(b131, ctx.path("gse131013_beta.parquet"))
    ctx.dump("gse131013_labels.pkl", pd.DataFrame({"y": lab2, "group": grp2}))
    tissue = md2.get("tissue", pd.Series(dtype=str)).value_counts().to_dict()
    info["GSE131013"] = {"n_samples": int(b131.shape[1]), "groups": grp2.value_counts().to_dict(),
                         "series_matrix_tissue": tissue, "n_series_samples": int(len(md2)),
                         "beta_range": [float(np.nanmin(b131.values)), float(np.nanmax(b131.values))],
                         "missing_network_prioritized": sorted(set(_comparator(ctx)) - set(b131.index))}
    for k in ("GSE48684", "GSE131013"):
        lo, hi = info[k]["beta_range"]
        assert -0.01 <= lo and hi <= 1.01, f"{k} values are not beta values: {lo}..{hi}"
    return info


def stage_classify(ctx: Ctx) -> dict:
    cfg = ctx.cfg
    c = cfg["classify"]
    prim = cfg["designs"]["primary"]
    seed = cfg["seed"]
    coad = ctx.load_info("coad")
    beta = D.load_frame(ctx.path("coad_beta.parquet"))      # discovery (imputed on all samples)
    samples = coad["tumors"] + coad["normals"]
    raw = _model_matrix(ctx, samples)                         # evaluation (imputed in-fold)
    groups = _patients(samples)
    y = np.r_[np.ones(len(coad["tumors"])), np.zeros(len(coad["normals"]))].astype(int)
    sept9 = [p for p in ctx.load("sept9.pkl") if p in beta.index]
    # selection funnel, last two steps: blood compatibility, then the panel rule
    spec = _specific(ctx, prim)
    bref = _blood_reference(ctx, _promoter_probes(ctx, raw.index))
    ok = B.blood_compatible(bref.reindex(spec.index), beta.loc[spec.index, coad["tumors"]].median(axis=1), cfg)
    panels = {"blood_compatible_11": C.delta_panel(spec[ok.values], c["panel_size"]),
              "screened_11": C.delta_panel(spec, c["panel_size"]),
              "network_prioritized_11": [p for p in _comparator(ctx) if p in beta.index],
              "SEPT9_promoter": sept9}
    geo = {"GSE48684": (D.load_frame(ctx.path("gse48684_beta.parquet")), ctx.load("gse48684_labels.pkl")),
           "GSE131013": (D.load_frame(ctx.path("gse131013_beta.parquet")), ctx.load("gse131013_labels.pkl"))}

    def _locked():
        rows, ext_rows, scores = [], [], {}
        for pname, panel in panels.items():
            X = raw.loc[panel, samples].T
            for model in c["models"]:
                cvdf, oof = C.locked_panel_cv(X, y, model, cfg, seed, groups=groups)
                m = C.metrics(y, np.nanmean(oof, axis=0), n_boot=c["n_boot_auc"], rs=D.rng(cfg, 11))
                rows.append({"panel": pname, "model": model, "eval": f"TCGA {c['cv_repeats']}x{c['cv_folds']}-fold CV (locked panel)",
                             **{k: m[k] for k in ("n", "tp", "fn", "tn", "fp", "sensitivity", "specificity", "accuracy", "auc")},
                             "auc_ci": m["auc_ci"], "sensitivity_ci": m["sensitivity_ci"],
                             "specificity_ci": m["specificity_ci"],
                             "acc_sd_over_repeats": float(cvdf["accuracy"].std()),
                             "auc_sd_over_repeats": float(cvdf["auc"].std())})
                scores[(pname, model, "TCGA_oof")] = (y, np.nanmean(oof, axis=0))
                sp = C.single_split(X, y, model, cfg, seed, D.rng(cfg, 12), groups=groups)
                rows.append({"panel": pname, "model": model, "eval": "TCGA single 70/30 split",
                             **{k: sp[k] for k in ("n", "tp", "fn", "tn", "fp", "sensitivity", "specificity", "accuracy", "auc")},
                             "auc_ci": sp["auc_ci"], "sensitivity_ci": sp["sensitivity_ci"],
                             "specificity_ci": sp["specificity_ci"]})
                exts = {}
                for gname, (gb, gl) in geo.items():
                    miss = [p for p in panel if p not in gb.index]
                    if miss:
                        ext_rows.append({"panel": pname, "model": model, "dataset": gname,
                                         "subset": "all", "note": f"missing CpGs: {miss}"})
                        continue
                    Xg = gb.loc[panel].T
                    exts[gname] = (Xg, gl.loc[Xg.index, "y"].values)
                    if gname == "GSE131013":
                        for sub, keep in {"tumor_vs_adjacent": gl["group"] != "healthy_mucosa",
                                          "tumor_vs_healthy": gl["group"] != "adjacent_normal"}.items():
                            idx = gl.index[keep.values]
                            exts[f"{gname}:{sub}"] = (Xg.loc[idx], gl.loc[idx, "y"].values)
                res, _ = C.frozen_external(X, y, exts, model, cfg, seed, D.rng(cfg, 13))
                for k, r in res.items():
                    gname, _, sub = k.partition(":")
                    ext_rows.append({"panel": pname, "model": model, "dataset": gname, "subset": sub or "all",
                                     **{kk: r.get(kk) for kk in ("n", "n_pos", "n_neg", "tp", "fn", "tn", "fp",
                                                                 "sensitivity", "specificity", "accuracy", "auc",
                                                                 "auc_ci", "sensitivity_ci", "specificity_ci",
                                                                 "unlabelled_n", "unlabelled_called_positive",
                                                                 "unlabelled_called_positive_ci")}})
                    if not sub:
                        lab = ~np.isnan(exts[k][1])
                        scores[(pname, model, gname)] = (exts[k][1][lab], r["scores"][lab])
            # training-free mean-beta score, threshold fixed on TCGA (Youden)
            sc = C.mean_beta_score(beta.loc[panel, samples], panel).values
            thr = C.youden_threshold(y, sc)
            m = C.metrics(y, sc, thr, n_boot=c["n_boot_auc"], rs=D.rng(cfg, 14))
            rows.append({"panel": pname, "model": "mean_beta", "eval": "TCGA (no training; threshold = Youden)",
                         **{k: m[k] for k in ("n", "tp", "fn", "tn", "fp", "sensitivity", "specificity", "accuracy", "auc")},
                         "auc_ci": m["auc_ci"], "threshold": thr})
            for gname, (gb, gl) in geo.items():
                if all(p in gb.index for p in panel):
                    s = C.mean_beta_score(gb.loc[panel], panel)
                    yy = gl.loc[s.index, "y"].values
                    lab = ~np.isnan(yy)
                    m = C.metrics(yy[lab], s.values[lab], thr, n_boot=c["n_boot_auc"], rs=D.rng(cfg, 15))
                    ext_rows.append({"panel": pname, "model": "mean_beta", "dataset": gname, "subset": "all",
                                     **{k: m[k] for k in ("n", "n_pos", "n_neg", "tp", "fn", "tn", "fp", "sensitivity",
                                                          "specificity", "accuracy", "auc", "auc_ci")}})
        return rows, ext_rows, scores

    rows, ext_rows, scores = ctx.cached("classify", "classify_locked.pkl", _locked)
    cv = pd.DataFrame(rows)
    ext = pd.DataFrame(ext_rows)
    ctx.table(cv, "classification_tcga.csv", index=False)
    ctx.table(ext, "classification_external.csv", index=False)
    ctx.dump("roc_scores.pkl", scores)

    log("  null panels")
    dmc = D.load_frame(ctx.path(f"dmc_{prim}.parquet"))
    pool = dmc.index[dmc["status"] == "hyper"]
    null = ctx.cached("classify", "classify_null.pkl", lambda: C.null_panels(
        raw, y, pool, len(panels["network_prioritized_11"]), c["n_null_panels"], D.rng(cfg, 21), seed,
        groups=groups))
    obs = {p: C.quick_cv_auc(raw.loc[panels[p], samples].T.to_numpy(), y, seed, groups=groups)
           for p in panels}
    pd.DataFrame({"null_auc": null}).to_csv(ctx.res / "tables" / "null_panel_auc.csv", index=False)
    null_info = {"pool": "hyper DMCs (design %s)" % prim, "n_pool": int(len(pool)),
                 "null_auc_median": float(np.median(null)),
                 "null_auc_q05_q95": np.percentile(null, [5, 95]).tolist(),
                 "observed": {p: {"auc": a, "empirical_p": float((1 + (null >= a).sum()) / (1 + len(null)))}
                              for p, a in obs.items()}}

    log("  nested CV (selection redone inside each outer fold)")
    nested = _nested(ctx, raw, y, groups, coad, bref)
    # determinism: rerun one locked-CV evaluation and compare hashes
    X = raw.loc[panels["blood_compatible_11"], samples].T
    h = [hashlib.md5(C.locked_panel_cv(X, y, "logreg", cfg, seed, groups=groups)[1].tobytes()).hexdigest()
         for _ in range(2)]
    return {"panels": panels, "null": null_info, "nested": nested,
            "n_specific": int(len(spec)), "n_specific_blood_compatible": int(ok.sum()),
            "determinism_hashes": h, "deterministic": h[0] == h[1]}


def _selection_half(ctx: Ctx, wb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Whole-blood donors split once: the first half defines blood compatibility, the second is
    simulated, so no blood sample both selects a CpG and scores it."""
    order = D.rng(ctx.cfg, 41).permutation(wb.columns)
    half = len(order) // 2
    return wb[order[:half]], wb[order[half:]]


def _promoter_probes(ctx: Ctx, probes) -> set:
    """Probes in the promoter definition of the panels: every CpG a fold could select."""
    ann = ctx.load("annot.pkl")
    groups = ctx.cfg["promoter_groups"][ctx.cfg["classify"]["panel_promoter_groups"]]
    g = A.genes_in(ann.loc[ann.index.intersection(probes)], groups)
    return set(g.index[g.map(len) > 0])


def _blood_reference(ctx: Ctx, probes) -> pd.DataFrame:
    """Blood methylation (B.blood_profile) of every probe a fold could select, from the
    selection half of GSE40279 and from GSE35069. Blood data take no part in the TCGA folds,
    so one profile serves every fold; only the tumour median is recomputed per fold."""
    cfg = ctx.cfg

    def _read():
        wb_path = D.download(cfg, cfg["urls"]["gse40279_matrix"])
        wb = D.series_table(wb_path, probes)
        c_path = D.download(cfg, cfg["urls"]["gse35069_matrix"])
        md_c, _ = D.parse_series_matrix(c_path, with_table=False)
        cells = D.series_table(c_path, probes)       # raw files stay for stage_blood
        wb_sel, _ = _selection_half(ctx, wb)
        return B.blood_profile(sorted(probes), wb_sel, cells, md_c["tissue/cell type"])

    # own checkpoint name: forcing classify must not re-download 656 blood profiles,
    # but a profile built for fewer probes than a fold can now select is rebuilt
    if ctx.path("blood_reference.pkl").exists():
        ref = ctx.load("blood_reference.pkl")
        if set(probes) <= set(ref.index):
            return ref
    ref = _read()
    ctx.dump("blood_reference.pkl", ref)
    return ref


def _nested(ctx: Ctx, beta: pd.DataFrame, y: np.ndarray, patients: np.ndarray, coad: dict,
            bref: pd.DataFrame) -> dict:
    """Outer 5-fold CV, grouped by patient, around the *entire* marker-discovery funnel
    (design B): missingness filter and imputation, DMCs, DEGs, integration, specificity
    screen, blood compatibility and the panel rule, all refitted on the training fold.
    The same folds also score the funnel without its blood step. beta is not imputed."""
    cfg = ctx.cfg
    c = cfg["classify"]
    pdef = c["panel_promoter_groups"]
    groups = cfg["promoter_groups"][pdef]
    ann = ctx.load("annot.pkl")
    counts = D.load_frame(ctx.path("rna_counts.parquet"))
    rna = ctx.load_info("rna")
    cond = pd.Series(["tumor"] * len(rna["tumors"]) + ["normal"] * len(rna["normals"]),
                     index=rna["tumors"] + rna["normals"])
    pm = ctx.load("probemap.pkl")
    patient = D.patient_of
    pan ={cc: D.load_frame(ctx.path("pancan") / f"{cc}.stats.parquet") for cc in cfg["screen_cancers"]}
    meta = {cc: m for cc, m in ctx.load_info("pancan")["cancers"].items() if cc in cfg["screen_cancers"]}
    comparator = set(_comparator(ctx))
    trace, memo = [], {}

    def select(train_cols, variant):
        key = frozenset(train_cols)          # logreg and mlp share the outer folds
        if key not in memo:                  # and a restarted run reuses finished folds
            # the tag keeps folds of the earlier sample-level, globally imputed funnel out
            h = hashlib.md5("|".join(["grouped-infold-v3-pydeseq2", *sorted(train_cols)]).encode()).hexdigest()[:12]
            memo[key] = ctx.cached("classify", f"nested_fold_funnel_{h}.pkl", lambda: _select(train_cols))
        trace.append(memo[key][1])
        return memo[key][0][variant]

    def _select(train_cols):
        train_set = set(train_cols)
        test_patients = {patient(s) for s in beta.columns if s not in train_set}
        tumors = [s for s in coad["tumors"] if s in train_set]
        normals = [s for s in coad["normals"] if s in train_set]
        assert not {patient(s) for s in tumors + normals} & test_patients
        b = M.filter_impute(beta[tumors + normals], cfg)      # fitted on the training fold only
        dmc = M.dmc_table(b, tumors, normals, cfg)
        # same procedure as stage_rna, on the training patients only: low-count filter, then PyDESeq2
        rna_cols = [s for s in cond.index if patient(s) not in test_patients]
        assert {"tumor", "normal"} <= set(cond.loc[rna_cols])
        filt = G.filter_genes(counts[rna_cols], cfg)
        log(f"    fold DESeq2 on {filt.shape[0]} genes, "
            f"{int((cond.loc[rna_cols] == 'tumor').sum())} T / {int((cond.loc[rna_cols] == 'normal').sum())} N")
        de = G.run_deseq2(filt, cond, cfg)
        de = G.to_symbols(de, pm)
        st = G.label_deg(de, cfg["deg"]["padj"], cfg["deg"]["min_abs_log2fc"])
        dmg = I.promoter_dmg(dmc, ann, groups)
        q = I.quadrants(dmg, st)
        cands = I.candidate_cpgs(dmc, ann, q["hyper_down"], groups)
        scr = S.screen(cands.index, pan, meta, cfg)
        spec = cands.loc[scr.index[scr["specific"]]]
        rec = {"n_dmc": int((dmc["status"] != "ns").sum()),
               "n_rna_samples": len(rna_cols), "n_genes_tested": int(filt.shape[0]),
               "n_deg_up": int((st == "up").sum()), "n_deg_down": int((st == "down").sum()),
               "n_hyper_down": len(q["hyper_down"]),
               "n_candidates": int(len(cands)), "n_specific": int(len(spec))}
        if spec.empty:
            return {"blood_compatible": [], "screened": []}, {**rec, "n_blood_compatible": 0, "panel": [],
                                                            "panel_screened": [], "network_prioritized_overlap": 0,
                                                            "network_prioritized_overlap_screened": 0}
        # training tumours only: the fold never sees the median of its test tumours
        ok = B.blood_compatible(bref.reindex(spec.index), b.loc[spec.index, tumors].median(axis=1), cfg)
        panels = {"blood_compatible": C.delta_panel(spec[ok.values], c["panel_size"]) if ok.any() else [],
                  "screened": C.delta_panel(spec, c["panel_size"])}
        rec.update({"n_blood_compatible": int(ok.sum()), "panel": panels["blood_compatible"],
                    "panel_screened": panels["screened"],
                    "network_prioritized_overlap": len(set(panels["blood_compatible"]) & comparator),
                    "network_prioritized_overlap_screened": len(set(panels["screened"]) & comparator)})
        log(f"    fold: {rec['n_dmc']} DMC, {len(spec)} specific, {rec['n_blood_compatible']} blood-compatible")
        return panels, rec

    out = {}
    for variant, slot in (("blood_compatible", out), ("screened", out.setdefault("screened", {}))):
        for model in ("logreg", "mlp"):
            trace.clear()
            res, _, _ = C.nested_cv(beta, y, lambda cols: select(cols, variant), model, cfg, cfg["seed"],
                                    D.rng(cfg, 31), groups=patients)
            slot[model] = {k: res[k] for k in ("n_scored", "tp", "fn", "tn", "fp", "sensitivity",
                                               "specificity", "accuracy", "auc", "auc_ci",
                                               "sensitivity_ci", "specificity_ci")}
            slot[model]["folds"] = list(trace)
            if model == "logreg":
                col = "panel" if variant == "blood_compatible" else "panel_screened"
                freq = pd.Series([p for t in trace for p in t[col]]).value_counts()
                slot["cpg_selection_frequency"] = freq.to_dict()
                if variant == "blood_compatible":
                    pd.DataFrame(trace).to_csv(ctx.res / "tables" / "nested_cv_folds.csv", index=False)
    return out


def stage_blood(ctx: Ctx) -> dict:
    """Leukocyte background of the panels, in-silico cfDNA dilution, real plasma cfDNA (EPIC)."""
    cfg = ctx.cfg
    b = cfg["blood"]
    c = cfg["classify"]
    prim = cfg["designs"]["primary"]
    seed = cfg["seed"]
    coad = ctx.load_info("coad")
    beta = D.load_frame(ctx.path("coad_beta.parquet"))
    spec = _specific(ctx, prim)
    panels = dict(ctx.load_info("classify")["panels"])
    probes = set(spec.index) | {p for v in panels.values() for p in v}

    def _blood():
        out = []
        for key in ("gse40279_matrix", "gse35069_matrix", "gse122126_epic_matrix"):
            path = D.download(cfg, cfg["urls"][key])
            md, _ = D.parse_series_matrix(path, with_table=False)
            out += [D.series_table(path, probes), md]
            keep_raw(ctx, path)
        wb, md_wb, cells, md_c, epic, md_m = out
        return wb, md_wb, cells, md_c["tissue/cell type"], epic, md_m

    wb, md_wb, cells, cell_type, epic, md_m = ctx.cached("blood", "blood_inputs.pkl", _blood)
    assert (md_wb["tissue"].str.lower() == "whole blood").all()
    wb_sel, wb_sim = _selection_half(ctx, wb)
    half = wb_sel.shape[1]
    sim_cases, sim_ctrl = wb_sim.iloc[:, : wb_sim.shape[1] // 2], wb_sim.iloc[:, wb_sim.shape[1] // 2:]

    keep = sorted(probes & set(beta.index))
    tum, nor = beta.loc[keep, coad["tumors"]], beta.loc[keep, coad["normals"]]
    for name, df in {"GSE40279": wb, "GSE35069": cells, "GSE122126": epic}.items():
        v = df.to_numpy()
        assert np.nanmin(v) >= 0 and np.nanmax(v) <= 1, f"{name}: values outside [0, 1]"
    prof = B.blood_profile(keep, wb_sel, cells, cell_type)
    prof.insert(0, "tcga_tumor_median", tum.median(axis=1).reindex(prof.index))
    prof.insert(1, "tcga_normal_median", nor.median(axis=1).reindex(prof.index))
    prof["wb_p95_holdout"] = wb_sim.quantile(0.95, axis=1).reindex(prof.index)
    prof["margin"] = prof["tcga_tumor_median"] - prof["blood_beta"]
    prof["blood_compatible"] = B.blood_compatible(prof, prof["tcga_tumor_median"], cfg)
    if ctx.path("blood_reference.pkl").exists():     # the nested folds read the same donors
        ref = ctx.load("blood_reference.pkl")["blood_beta"]
        both = prof.index.intersection(ref.index)
        assert np.allclose(prof.loc[both, "blood_beta"], ref[both], equal_nan=True), "blood reference drift"

    groups = B.moss_groups(md_m)
    for g in ("healthy_cfdna", "colon_cancer_cfdna", "colon_epithelium", "leukocytes"):
        cols = groups.index[(groups["group"] == g) & (groups["replicate_of"] == "")]
        prof[f"moss_{g}_median"] = epic.reindex(prof.index)[cols].median(axis=1)
    prof["on_epic"] = prof.index.isin(epic.index)

    ok = prof["blood_compatible"].reindex(spec.index).fillna(False).astype(bool)
    panels["network_prioritized_blood_compatible"] = [p for p in panels["network_prioritized_11"] if prof.loc[p, "blood_compatible"]]
    assert C.delta_panel(spec[ok], c["panel_size"]) == panels["blood_compatible_11"], "panel drift since classify"
    panels = {k: v for k, v in panels.items() if v}
    member = pd.DataFrame({k: prof.index.isin(v) for k, v in panels.items()}, index=prof.index)
    prof = prof.join(member.add_prefix("in_"))
    # CpGs outside the specific set (e.g. MAL, SUSD5) still get their annotated genes
    annotated = ctx.load("annot.pkl")["gene_group"].reindex(prof.index).map(
        lambda gg: ";".join(sorted({g for g, _ in gg})) or np.nan if isinstance(gg, list) else np.nan)
    prof["genes"] = spec["genes"].reindex(prof.index).fillna(annotated)
    ctx.table(prof[member.any(axis=1)], "blood_background_panels.csv")
    ctx.table(prof.loc[prof.index.isin(spec.index)], f"blood_background_specific_{prim}.csv")

    log("  in-silico dilution")
    sim = B.simulate_dilution(tum, sim_cases, sim_ctrl, panels, b["tumor_fractions"], b["read_depths"],
                              b["n_mixtures"], b["specificity"], D.rng(cfg, 42))
    ctx.table(sim, "blood_dilution.csv", index=False)
    lod = B.detection_limit(sim, b["lod_sensitivity"])
    ctx.table(lod, "blood_detection_limit.csv", index=False)

    log("  real cfDNA (GSE122126, EPIC)")
    cf = B.cfdna_scores(epic, groups, panels)
    ctx.table(cf, "cfdna_scores.csv", index=False)
    ctx.table(groups, "cfdna_groups.csv")

    log("  blood-compatible panel: TCGA CV and frozen external validation")
    samples = coad["tumors"] + coad["normals"]
    raw, patients = _model_matrix(ctx, samples), _patients(samples)
    y = np.r_[np.ones(len(coad["tumors"])), np.zeros(len(coad["normals"]))].astype(int)
    geo = {"GSE48684": (D.load_frame(ctx.path("gse48684_beta.parquet")), ctx.load("gse48684_labels.pkl")),
           "GSE131013": (D.load_frame(ctx.path("gse131013_beta.parquet")), ctx.load("gse131013_labels.pkl"))}
    rows = []
    for pname in ("network_prioritized_blood_compatible", "blood_compatible_11"):
        panel = panels.get(pname, [])
        if not panel:
            continue
        X = raw.loc[panel, samples].T
        for model in ("logreg", "mlp"):
            _, oof = C.locked_panel_cv(X, y, model, cfg, seed, groups=patients)
            m = C.metrics(y, np.nanmean(oof, axis=0), n_boot=c["n_boot_auc"], rs=D.rng(cfg, 43))
            rows.append({"panel": pname, "model": model, "dataset": "TCGA CV (locked panel)",
                         **{k: m.get(k) for k in _METRIC_KEYS}})
            exts = {g: (gb.loc[panel].T, gl.loc[gb.columns, "y"].values)
                    for g, (gb, gl) in geo.items() if all(p in gb.index for p in panel)}
            res, _ = C.frozen_external(X, y, exts, model, cfg, seed, D.rng(cfg, 44))
            for g, r in res.items():
                rows.append({"panel": pname, "model": model, "dataset": f"{g} (frozen TCGA model)",
                             **{k: r.get(k) for k in _METRIC_KEYS + ["unlabelled_n", "unlabelled_called_positive"]}})
    cls = pd.DataFrame(rows)
    ctx.table(cls, "blood_panel_classification.csv", index=False)

    pub = prof.loc[[p for p in _comparator(ctx) if p in prof.index]]
    colon = cf[(cf["group"] == "colon_cancer_cfdna") & (cf["replicate_of"] == "")]
    spike = cf[cf["group"] == "cfdna_spike_in"]
    return {"n_whole_blood": int(wb.shape[1]), "n_whole_blood_selection": int(half),
            "n_sorted_cells": int(cells.shape[1]), "cell_types": sorted(cell_type.unique()),
            "n_epic_samples": int(epic.shape[1]), "epic_groups": groups["group"].value_counts().to_dict(),
            "n_specific": int(len(spec)), "n_specific_blood_compatible": int(ok.sum()),
            "network_prioritized_blood_beta": pub["blood_beta"].round(3).to_dict(),
            "network_prioritized_blood_compatible": panels.get("network_prioritized_blood_compatible", []),
            "network_prioritized_not_on_epic": sorted(set(_comparator(ctx)) - set(epic.index)),
            "panels": panels,
            "lowest_fraction_detected": {f"{r.panel}|depth{r.depth}": r.lowest_fraction_detected
                                         for r in lod.itertuples()},
            "colon_cfdna_z": {f"{r.panel}|{r.title}": round(r.z_vs_healthy, 2) for r in colon.itertuples()},
            "spike_in_z": {f"{r.panel}|{r.title}|{r.spike_fraction}": round(r.z_vs_healthy, 2)
                           for r in spike.itertuples()},
            "classification": cls.drop(columns=[c_ for c_ in cls.columns if c_.endswith("_ci")]).to_dict("records")}


_METRIC_KEYS = ["n", "n_pos", "n_neg", "tp", "fn", "tn", "fp", "sensitivity", "specificity",
                "accuracy", "auc", "auc_ci", "sensitivity_ci", "specificity_ci"]


def stage_report(ctx: Ctx) -> dict:
    return R.build(ctx)


# --------------------------------------------------------------------------- driver
def main(argv=None) -> None:
    import warnings
    from sklearn.exceptions import ConvergenceWarning
    warnings.filterwarnings("once", category=ConvergenceWarning)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(Path(__file__).with_name("config.yaml")))
    ap.add_argument("--workdir")
    ap.add_argument("--resultsdir")
    ap.add_argument("--stage", choices=STAGES, help="run only this stage")
    ap.add_argument("--from", dest="start", choices=STAGES, help="run from this stage on")
    ap.add_argument("--until", choices=STAGES, help="stop after this stage")
    ap.add_argument("--force", action="store_true", help="recompute selected stages")
    args = ap.parse_args(argv)
    cfg = D.load_config(args.config, args.workdir, args.resultsdir)
    if args.stage:
        todo = [args.stage]
    else:
        i = STAGES.index(args.start) if args.start else 0
        j = STAGES.index(args.until) + 1 if args.until else len(STAGES)
        todo = STAGES[i:j]
    ctx = Ctx(cfg, set(todo) if args.force else set())
    shutil.copy(args.config, ctx.res / "config_used.yaml")
    D.write_json(D.env_versions(), ctx.res / "environment.json")
    import subprocess
    freeze = subprocess.run([sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True)
    (ctx.res / "logs" / "pip_freeze.txt").write_text(freeze.stdout)
    for stage in todo:
        if ctx.done(stage):
            log(f"{stage}: cached")
            continue
        log(f"{stage}: start")
        t0 = time.time()
        info = globals()[f"stage_{stage}"](ctx)
        info["seconds"] = round(time.time() - t0, 1)
        ctx.mark(stage, info)
        log(f"{stage}: done in {info['seconds']}s")


if __name__ == "__main__":
    main()
