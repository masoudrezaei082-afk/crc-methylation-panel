#!/usr/bin/env python
"""End-to-end smoke run of run_all.py on small synthetic inputs (no internet access).

Fabricates every raw file under <tmp>/work/raw with the names download() expects,
fabricates the blood and plasma series matrices, replaces the enrichment web service (Enrichr) with a stub and runs all stages. It checks that the plumbing works, not the biology.

    python tests/smoke_pipeline.py [outdir]
"""

import gzip
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import run_all  # noqa: E402
from crcmeth import integrate as I  # noqa: E402

RS = np.random.default_rng(0)
PANEL = {  # probe -> (gene, region, island) as in GPL13534
    "cg05142617": ("MAL", "TSS1500", "N_Shore"), "cg01217071": ("AOX1", "TSS1500", "N_Shore"),
    "cg00329270": ("KY", "TSS200", "Island"), "cg03158581": ("HAND1", "TSS1500", "Island"),
    "cg16683920": ("TDRD10", "TSS1500", "Island"), "cg10775792": ("PCP4L1", "TSS200", "N_Shore"),
    "cg27223727": ("PCP4L1", "TSS200", "Island"), "cg15496063": ("CYS1", "TSS1500", "Island"),
    "cg08385249": ("LIMS2", "TSS1500", "Island"), "cg22657780": ("ZNF385B", "1stExon", "Island"),
    "cg04545136": ("SUSD5", "TSS200", "Island")}


def barcodes(n_t, n_n, prefix):
    t = [f"TCGA-{prefix}-{i:04d}-01A" for i in range(n_t)]
    n = [f"TCGA-{prefix}-{i:04d}-11A" for i in range(n_n)]
    return t, n


def write_tsv_gz(df, path, index_label):
    with gzip.open(path, "wt") as fh:
        df.to_csv(fh, sep="\t", index_label=index_label, float_format="%.4f", na_rep="NA")


def fabricate(raw: Path, n_probes=3000):
    raw.mkdir(parents=True, exist_ok=True)
    genes = [f"GENE{i}" for i in range(300)] + sorted({g for g, _, _ in PANEL.values()}) + ["SEPT9"]
    probes = list(PANEL) + [f"cg9{i:07d}" for i in range(n_probes - len(PANEL))]
    regions = ["TSS1500", "TSS200", "Body", "5'UTR", "1stExon", "3'UTR"]
    islands = ["Island", "N_Shore", "S_Shore", "OpenSea", ""]
    rows = []
    for p in probes:
        if p in PANEL:
            g, r, isl = PANEL[p]
        else:
            g, r, isl = RS.choice(genes), RS.choice(regions), RS.choice(islands)
        chrom = "X" if RS.random() < 0.02 else str(RS.integers(1, 23))
        rows.append((p, chrom, int(RS.integers(1, 1e8)), g, r, isl))
    ann = pd.DataFrame(rows, columns=["IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name",
                                      "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"])
    with gzip.open(raw / "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz", "wt") as fh:
        fh.write("Illumina, Inc.\n[Heading]\n[Assay]\n")
        ann.to_csv(fh, index=False)
    masked = set(RS.choice(probes[len(PANEL):], 50, replace=False))
    with gzip.open(raw / "HM450.hg38.mask.tsv.gz", "wt") as fh:
        pd.DataFrame({"Probe_ID": probes, "mask": "", "maskUniq": "",
                      "M_general": [str(p in masked).upper() for p in probes]}).to_csv(fh, sep="\t", index=False)
    pm = pd.DataFrame({"id": [f"ENSG{i:011d}.1" for i in range(len(genes))], "gene": genes,
                       "chrom": "chr1", "chromStart": RS.integers(1e6, 2e8, len(genes)),
                       "chromEnd": 0, "strand": RS.choice(["+", "-"], len(genes))})
    pm["chromEnd"] = pm["chromStart"] + 5000
    pm.to_csv(raw / "gencode_v36_probemap.tsv", sep="\t", index=False)

    # methylation: hyper probes in tumours on a set of "silenced" genes (incl. the comparator panel)
    silenced = set(genes[:60]) | {g for g, _, _ in PANEL.values()}
    hyper = np.array([ann.UCSC_RefGene_Name[i] in silenced and ann.UCSC_RefGene_Group[i] in
                      ("TSS1500", "TSS200", "1stExon") for i in range(len(probes))])
    base = RS.beta(2, 6, len(probes))

    def meth(n_t, n_n, shift_mask, shift):
        t = np.clip(base[:, None] + shift_mask[:, None] * shift + RS.normal(0, 0.05, (len(probes), n_t)), 0, 1)
        n = np.clip(base[:, None] + RS.normal(0, 0.05, (len(probes), n_n)), 0, 1)
        m = np.c_[t, n]
        m[RS.random(m.shape) < 0.003] = np.nan
        return m

    t, n = barcodes(70, 12, "AA")
    write_tsv_gz(pd.DataFrame(meth(70, 12, hyper, 0.4), index=probes, columns=t + n),
                 raw / "TCGA-COAD.methylation450.tsv.gz", "Composite Element REF")
    shared = hyper & (RS.random(len(probes)) < 0.5)        # half of COAD hyper probes pan-cancer
    for k, c in enumerate(["PAAD", "BRCA", "BLCA", "LIHC", "THCA", "KIRC", "LUSC", "SARC",
                           "ESCA", "CESC", "LUAD"]):
        other = shared
        t2, n2 = barcodes(15, 12, f"{k:02d}")
        write_tsv_gz(pd.DataFrame(meth(15, 12, other, 0.4), index=probes, columns=t2 + n2),
                     raw / f"TCGA-{c}.methylation450.tsv.gz", "Composite Element REF")

    # RNA: silenced genes down in tumours
    rt, rn = barcodes(60, 12, "AA")
    mu = RS.lognormal(5, 1, len(genes))
    fold = np.where(np.isin(genes, list(silenced)), 0.15, 1.0)
    cnt = np.c_[RS.poisson(mu[:, None] * fold[:, None], (len(genes), 60)),
                RS.poisson(mu[:, None], (len(genes), 12))]
    write_tsv_gz(pd.DataFrame(np.log2(cnt + 1), index=pm["id"], columns=rt + rn),
                 raw / "TCGA-COAD.star_counts.tsv.gz", "Ensembl_ID")

    # GEO
    gsm = [f"GSM{i}" for i in range(40)]
    status = ["cancer"] * 15 + ["adenoma"] * 10 + ["normal-C"] * 8 + ["normal-H"] * 7
    m = meth(15, 25, hyper, 0.35)
    tab = pd.DataFrame(m, index=probes, columns=gsm)
    with gzip.open(raw / "GSE48684_series_matrix.txt.gz", "wt") as fh:
        fh.write("!Series_title\t\"x\"\n")
        fh.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
        fh.write("!Sample_characteristics_ch1\t" + "\t".join(f'"disease status: {s}"' for s in status) + "\n")
        fh.write("!series_matrix_table_begin\n")
        tab.to_csv(fh, sep="\t", index_label="ID_REF", float_format="%.4f", na_rep="null")
        fh.write("!series_matrix_table_end\n")
    codes = [f"A{i:04d}_T" for i in range(20)] + [f"A{i:04d}_N" for i in range(20)] + \
            [f"B{i:04d}_M" for i in range(10)]
    m2 = meth(20, 30, hyper, 0.35)
    cols = {}
    for j, c in enumerate(codes):
        cols[c] = m2[:, j]
        cols[f"{c}.Detection.Pval"] = 0.0
    write_tsv_gz(pd.DataFrame(cols, index=probes), raw / "GSE131013_normalized_matrix.txt.gz", "ID_REF")
    with gzip.open(raw / "GSE131013_series_matrix.txt.gz", "wt") as fh:
        g2 = [f"GSM9{i}" for i in range(len(codes))]
        fh.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in g2) + "\n")
        fh.write("!Sample_characteristics_ch1\t" + "\t".join(
            f'"tissue: {"Tumor" if c.endswith("T") else "Normal"}"' for c in codes) + "\n")
        fh.write("!series_matrix_table_begin\n\"ID_REF\"\t" + "\t".join(f'"{g}"' for g in g2) + "\n")
        fh.write("!series_matrix_table_end\n")

    # blood and plasma: half of the hyper probes are unmethylated in leukocytes
    low = hyper & (RS.random(len(probes)) < 0.5)
    blood = lambda n: np.clip(np.where(low, 0.02, base)[:, None] + RS.normal(0, 0.01, (len(probes), n)), 0, 1)
    series_matrix(raw / "GSE40279_series_matrix.txt.gz", probes, blood(40), {"tissue": ["whole blood"] * 40})
    cells = ["Whole blood", "PBMC", "Granulocytes", "CD4+ T cells", "CD8+ T cells", "CD14+ Monocytes"] * 2
    series_matrix(raw / "GSE35069_series_matrix.txt.gz", probes, blood(12), {"tissue/cell type": cells})
    titles = [f"cfDNA_healthy_{i}" for i in range(6)] + [f"cfDNA_CRC_{i}" for i in range(3)] + \
             ["colon_epi_1", "leukocytes_1"]
    source = ["cfDNA, healthy"] * 6 + ["cfDNA, colon adenocarcinoma"] * 3 + ["Colon epithelial cells", "Leukocytes"]
    epic = np.c_[blood(6), 0.7 * blood(3) + 0.3 * meth(3, 0, hyper, 0.4), meth(1, 0, hyper, 0.4), blood(1)]
    series_matrix(raw / "GSE122126-GPL21145_series_matrix.txt.gz", probes, epic, {}, title=titles, source_name_ch1=source)


def series_matrix(path, probes, values, chars, **fields):
    """Minimal GEO series matrix: sample fields, `key: value` characteristics and the table."""
    gsm = [f"GSM{path.name[3:8]}{i:03d}" for i in range(values.shape[1])]
    with gzip.open(path, "wt") as fh:
        fh.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
        for k, v in fields.items():
            fh.write(f"!Sample_{k}\t" + "\t".join(f'"{x}"' for x in v) + "\n")
        for k, v in chars.items():
            fh.write("!Sample_characteristics_ch1\t" + "\t".join(f'"{k}: {x}"' for x in v) + "\n")
        fh.write("!series_matrix_table_begin\n")
        pd.DataFrame(values, index=probes, columns=gsm).to_csv(fh, sep="\t", index_label="ID_REF",
                                                               float_format="%.4f", na_rep="null")
        fh.write("!series_matrix_table_end\n")


def stub_services():
    def enrichment(*a, **k):
        raise RuntimeError("enrichment stubbed in smoke test")

    I.enrichment = enrichment

    # STRING, UCSC and JASPAR are web services: replace them by seeded random edges
    def string_edges(genes, cfg, cache_dir, required_score=None):
        g = sorted(set(genes))
        rs = np.random.default_rng(len(g) + (required_score or 0))
        pairs = {tuple(sorted(rs.choice(g, 2, replace=False))) for _ in range(3 * len(g))}
        return pd.DataFrame([(a, b, float(rs.uniform(0.4, 1))) for a, b in sorted(pairs)],
                            columns=["a", "b", "weight"])

    def tf_edges(genes, probemap, cfg, cache_dir):
        g = sorted(set(genes))
        rs = np.random.default_rng(7)
        rows = [(f"TF{t}", x) for x in g for t in rs.choice(40, 3, replace=False)]
        return pd.DataFrame(rows, columns=["tf", "gene"]).assign(weight=1.0)

    def request(*a, **k):
        raise RuntimeError("web services stubbed in smoke test")

    run_all.N.string_edges, run_all.N.tf_edges, run_all.N._request = string_edges, tf_edges, request


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(tempfile.mkdtemp(prefix="crcmeth_smoke_"))
    fabricate(out / "work" / "raw")
    cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
    cfg["keep_raw"] = True
    cfg["design_a_per_source"] = 10
    cfg["specificity"]["min_normals_for_test"] = 5
    cfg["classify"].update(cv_repeats=2, cv_folds=5, n_null_panels=20, n_boot_auc=50,
                           rf_trees=50, mlp_max_iter=300)
    cfg["deg"]["n_cpus"] = 1
    cfg["network"].update(n_bootstrap=5, beta_grid_step=0.5)
    cfg_path = out / "config_smoke.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    stub_services()
    run_all.main(["--config", str(cfg_path), "--workdir", str(out / "work"),
                  "--resultsdir", str(out / "results")])
    # resume path: a crashed classify stage reuses its checkpoints
    (out / "work" / "cache" / "classify.done").unlink()
    assert list((out / "work" / "cache").glob("nested_fold_*.pkl"))
    run_all.main(["--config", str(cfg_path), "--workdir", str(out / "work"),
                  "--resultsdir", str(out / "results"), "--from", "classify"])
    res = out / "results"
    for f in ("SUMMARY.md", "metrics.json",
              "tables/traceability_network_prioritized_panel.csv", "tables/classification_external.csv"):
        assert (res / f).exists(), f
    print("SMOKE OK", res)


if __name__ == "__main__":
    main()
