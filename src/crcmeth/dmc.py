"""Probe QC, imputation and differential methylation (Wilcoxon rank-sum + BH)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def probe_qc(beta: pd.DataFrame, annot: pd.DataFrame, mask: set[str], cfg: dict,
             missing: bool = True) -> pd.DataFrame:
    """Annotation, sex-chromosome and mask filters, which depend on no sample; then, unless
    `missing` is False, the missingness filter, which depends on the samples given."""
    q = cfg["probe_qc"]
    keep = beta.index.isin(annot.index)
    if q["drop_sex_chrom"]:
        sex = set(annot.index[annot["chrom"].isin(["X", "Y"])])
        keep &= ~beta.index.isin(sex)
    if q["use_zhou_mask"]:
        keep &= ~beta.index.isin(mask)
    out = beta.loc[keep]
    return missing_filter(out, q["max_missing_frac"]) if missing else out


def missing_filter(beta: pd.DataFrame, max_frac: float) -> pd.DataFrame:
    return beta.loc[beta.isna().mean(axis=1) <= max_frac]


def filter_impute(beta: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Missingness filter and imputation, both fitted on exactly the columns given; inside
    cross-validation pass the training columns only."""
    q = cfg["probe_qc"]
    return impute_rows(missing_filter(beta, q["max_missing_frac"]), q["impute"]).astype(np.float32)


def impute_rows(beta: pd.DataFrame, method: str = "median") -> pd.DataFrame:
    """Row-wise fill from the columns given. Label-free, but the fill values depend on every
    column, so inside cross-validation it must see the training columns only."""
    fill = beta.median(axis=1) if method == "median" else beta.mean(axis=1)
    values = beta.to_numpy(copy=True)
    rows, cols = np.where(np.isnan(values))
    values[rows, cols] = fill.to_numpy()[rows]
    return pd.DataFrame(values, index=beta.index, columns=beta.columns)


def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaN p-values stay NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def dmc_table(beta: pd.DataFrame, tumor: list[str], normal: list[str], cfg: dict) -> pd.DataFrame:
    """Per-probe tumour-vs-normal statistics and DMC status (hyper/hypo/ns)."""
    d = cfg["dmc"]
    t = beta[tumor].to_numpy(np.float64)
    n = beta[normal].to_numpy(np.float64)
    p = np.empty(beta.shape[0])
    step = d["chunk_rows"]
    for i in range(0, beta.shape[0], step):
        p[i:i + step] = mannwhitneyu(t[i:i + step], n[i:i + step], axis=1,
                                     method="asymptotic", nan_policy="omit").pvalue
    out = pd.DataFrame({"mean_t": np.nanmean(t, axis=1), "mean_n": np.nanmean(n, axis=1),
                        "median_t": np.nanmedian(t, axis=1), "median_n": np.nanmedian(n, axis=1),
                        "p": p}, index=beta.index)
    out["delta"] = out["mean_t"] - out["mean_n"]
    out["fdr"] = bh(out["p"].to_numpy())
    sig = (out["fdr"] < d["fdr"]) & (out["delta"].abs() > d["min_abs_delta_beta"])
    out["status"] = np.where(sig & (out["delta"] > 0), "hyper",
                             np.where(sig, "hypo", "ns"))
    out.attrs.update(n_tumor=len(tumor), n_normal=len(normal))
    return out
