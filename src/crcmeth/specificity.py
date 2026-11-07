"""Pan-cancer streaming statistics and the CRC-specificity screen."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .data import NORMAL_CODE, TUMOR_CODE, header_columns, iter_matrix, select_samples
from .dmc import bh


def pancancer_stats(path, cfg: dict, keep_probes=frozenset(), n_extract: int = 0,
                    rng: np.random.Generator | None = None):
    """Stream one TCGA 450K matrix and return

    stats     per-probe tumour-vs-own-normal statistics (all probes)
    kept      tumour+normal values for `keep_probes` (figures, absolute rule)
    extracted all probes for `n_extract` random normals (Design A controls)
    """
    cols = header_columns(path)
    tumors = select_samples(cols, TUMOR_CODE)
    normals = select_samples(cols, NORMAL_CODE)
    extract = []
    if n_extract:
        extract = sorted(rng.choice(normals, size=min(n_extract, len(normals)),
                                    replace=False).tolist())
    max_na = cfg["probe_qc"]["max_missing_frac"]
    stats, kept, ext = [], [], []
    for chunk in iter_matrix(path, tumors + normals, cfg["dmc"]["chunk_rows"]):
        t = chunk[tumors].to_numpy(np.float64)
        n = chunk[normals].to_numpy(np.float64)
        bad = (np.isnan(t).mean(1) > max_na) | (np.isnan(n).mean(1) > max_na)
        p = np.full(len(chunk), np.nan)
        if len(normals) >= 2:
            ok = ~bad
            if ok.any():
                p[ok] = mannwhitneyu(t[ok], n[ok], axis=1, method="asymptotic",
                                     nan_policy="omit").pvalue
        with np.errstate(all="ignore"):
            stats.append(pd.DataFrame({
                "mean_t": np.nanmean(t, 1), "mean_n": np.nanmean(n, 1),
                "median_t": np.nanmedian(t, 1), "median_n": np.nanmedian(n, 1),
                "p": p}, index=chunk.index))
        sel = chunk.index.isin(keep_probes)
        if sel.any():
            kept.append(chunk.loc[sel])
        if extract:
            ext.append(chunk[extract])
    st = pd.concat(stats)
    st["delta"] = st["mean_t"] - st["mean_n"]
    st["fdr"] = bh(st["p"].to_numpy())
    st.attrs.update(n_tumor=len(tumors), n_normal=len(normals))
    kept_df = pd.concat(kept) if kept else pd.DataFrame(columns=tumors + normals)
    ext_df = pd.concat(ext) if ext else None
    meta = {"n_tumor": len(tumors), "n_normal": len(normals), "tumors": tumors,
            "normals": normals, "extracted_normals": extract}
    return st, kept_df, ext_df, meta


def screen(candidates, pancan: dict[str, pd.DataFrame], meta: dict[str, dict],
           cfg: dict) -> pd.DataFrame:
    """Rows = candidate CpGs, columns = per-cancer exclusion flags + 'specific'."""
    d, s = cfg["dmc"], cfg["specificity"]
    out = pd.DataFrame(index=pd.Index(sorted(candidates), name="probe"))
    for cancer, st in pancan.items():
        st = st.reindex(out.index)
        dmc = (st["fdr"] < d["fdr"]) & (st["delta"].abs() > d["min_abs_delta_beta"])
        flag = dmc.fillna(False)
        if s["rule"] == "dmc_or_high":
            flag |= (st["median_t"] > s["absolute_max_tumor_median"]).fillna(False)
        out[f"{cancer}_excluded"] = flag
        out[f"{cancer}_delta"] = st["delta"]
        out[f"{cancer}_median_t"] = st["median_t"]
    flags = [c for c in out.columns if c.endswith("_excluded")]
    out["n_cancers_excluding"] = out[flags].sum(axis=1)
    out["specific"] = out["n_cancers_excluding"] == 0
    out.attrs["underpowered"] = sorted(c for c, m in meta.items()
                                       if m["n_normal"] < s["min_normals_for_test"])
    return out
