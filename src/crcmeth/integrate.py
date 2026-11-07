"""Methylation x expression integration, genomic-region summaries, enrichment."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .annot import ISLAND_ORDER, REGION_ORDER, genes_in


def region_summary(dmc: pd.DataFrame, annot: pd.DataFrame, promoter_groups) -> dict:
    """Counts behind Fig 1 A-D (all / promoter / CGI / CGI promoter; island; gene region)."""
    a = annot.loc[dmc.index]
    sig = dmc["status"] != "ns"
    prom = genes_in(a, promoter_groups).map(len) > 0
    isl = a["island"] == "Island"
    out = {}
    for name, mask in {"all": sig, "promoter": sig & prom, "CGI": sig & isl,
                       "CGI_promoter": sig & prom & isl}.items():
        out[f"fig1A_{name}"] = dmc.loc[mask, "status"].value_counts().to_dict()
    out["fig1B_island"] = {
        k: dmc.loc[sig & (a["island"] == k), "status"].value_counts().to_dict()
        for k in ISLAND_ORDER}
    out["fig1C_region"] = {
        r: dmc.loc[sig & a["gene_group"].map(lambda gg: any(x == r for _, x in gg)),
                   "status"].value_counts().to_dict() for r in REGION_ORDER}
    genes = a["gene_group"].map(lambda gg: sorted({g for g, _ in gg}))
    out["fig1D_genes"] = {s: int(len(set().union(*genes[dmc["status"] == s])) if
                                 (dmc["status"] == s).any() else 0)
                          for s in ("hyper", "hypo")}
    return out


def promoter_dmg(dmc: pd.DataFrame, annot: pd.DataFrame, promoter_groups) -> pd.DataFrame:
    """Gene-level promoter methylation status from >=1 promoter DMC (stated rule)."""
    pg = genes_in(annot.loc[dmc.index], promoter_groups)
    rows = []
    for probe, genes in pg.items():
        st = dmc.at[probe, "status"]
        if st == "ns":
            continue
        rows.extend((g, probe, st) for g in genes)
    long = pd.DataFrame(rows, columns=["gene", "probe", "status"])
    if long.empty:
        return pd.DataFrame(columns=["n_hyper", "n_hypo", "status"])
    tab = long.pivot_table(index="gene", columns="status", values="probe",
                           aggfunc="count", fill_value=0)
    tab = tab.reindex(columns=["hyper", "hypo"], fill_value=0)
    tab.columns = ["n_hyper", "n_hypo"]
    tab["status"] = np.select([(tab.n_hyper > 0) & (tab.n_hypo == 0),
                               (tab.n_hypo > 0) & (tab.n_hyper == 0)],
                              ["hyper", "hypo"], "mixed")
    return tab


def quadrants(dmg: pd.DataFrame, deg_label: pd.Series) -> dict:
    out = {}
    for m in ("hyper", "hypo"):
        genes = set(dmg.index[dmg["status"] == m])
        out[f"{m}_genes"] = len(genes)
        out[f"{m}_genes_with_expression"] = len(genes & set(deg_label.index))
        for e in ("up", "down"):
            out[f"{m}_{e}"] = sorted(genes & set(deg_label.index[deg_label == e]))
    return out


def candidate_cpgs(dmc: pd.DataFrame, annot: pd.DataFrame, genes, promoter_groups) -> pd.DataFrame:
    """Hypermethylated promoter DMCs that annotate to one of `genes`."""
    genes = set(genes)
    pg = genes_in(annot.loc[dmc.index], promoter_groups)
    hit = pg.map(lambda gs: sorted(genes.intersection(gs)))
    keep = (dmc["status"] == "hyper") & (hit.map(len) > 0)
    out = dmc.loc[keep, ["delta", "fdr", "mean_t", "mean_n"]].copy()
    out["genes"] = hit[keep].map(";".join)
    return out


def matched_spearman(beta: pd.DataFrame, expr: pd.DataFrame, cpg_gene: pd.Series) -> pd.DataFrame:
    """Spearman rho between CpG beta and log expression over patients with both assays.

    beta / expr columns must already be patient ids; cpg_gene maps probe -> gene."""
    pts = sorted(set(beta.columns) & set(expr.columns))
    rows = []
    for probe, gene in cpg_gene.items():
        if probe not in beta.index or gene not in expr.index:
            continue
        rho, p = spearmanr(beta.loc[probe, pts], expr.loc[gene, pts], nan_policy="omit")
        rows.append((probe, gene, rho, p))
    return pd.DataFrame(rows, columns=["probe", "gene", "rho", "p"]).assign(n=len(pts))


def enrichment(genes, background, libraries, outdir=None) -> pd.DataFrame:
    """Hypergeometric over-representation against Enrichr libraries (DAVID substitute)."""
    import gseapy

    frames = []
    for lib in libraries:
        gmt = gseapy.get_library(name=lib, organism="Human")
        res = gseapy.enrich(gene_list=list(genes), gene_sets=gmt,
                            background=list(background), outdir=outdir, verbose=False)
        df = res.results.copy()
        df["library"] = lib
        frames.append(df)
    return pd.concat(frames, ignore_index=True)
