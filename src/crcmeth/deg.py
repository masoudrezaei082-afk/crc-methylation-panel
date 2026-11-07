"""RNA-seq differential expression: DESeq2 (pydeseq2), used in the main analysis and inside every nested-CV fold."""

from __future__ import annotations

import numpy as np
import pandas as pd


def counts_from_xena(log2_counts: pd.DataFrame) -> pd.DataFrame:
    """Xena stores log2(count + 1); recover integer counts (genes x samples)."""
    return np.rint(np.exp2(log2_counts.astype(np.float64)) - 1).clip(lower=0).astype(np.int64)


def filter_genes(counts: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    d = cfg["deg"]
    need = max(3, int(np.ceil(d["min_samples_frac"] * counts.shape[1])))
    return counts.loc[(counts >= d["min_count"]).sum(axis=1) >= need]


def run_deseq2(counts: pd.DataFrame, condition: pd.Series, cfg: dict) -> pd.DataFrame:
    """counts: genes x samples; condition: sample -> 'tumor' | 'normal'."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    inference = DefaultInference(n_cpus=cfg["deg"]["n_cpus"])
    meta = pd.DataFrame({"condition": condition.loc[counts.columns].values},
                        index=counts.columns)
    dds = DeseqDataSet(counts=counts.T, metadata=meta, design="~condition",
                       refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "tumor", "normal"],
                       inference=inference, quiet=True)
    stats.summary()
    res = stats.results_df.rename(columns={"log2FoldChange": "log2fc", "pvalue": "p"})
    return res[["baseMean", "log2fc", "lfcSE", "p", "padj"]]


def label_deg(res: pd.DataFrame, padj: float, min_abs_log2fc: float) -> pd.Series:
    sig = (res["padj"] < padj) & (res["log2fc"].abs() > min_abs_log2fc)
    return pd.Series(np.where(sig & (res["log2fc"] > 0), "up",
                              np.where(sig, "down", "ns")), index=res.index)


def to_symbols(res: pd.DataFrame, probemap: pd.DataFrame) -> pd.DataFrame:
    """Ensembl -> symbol; duplicated symbols keep the row with smallest padj."""
    out = res.join(probemap["gene"], how="inner")
    out = out.sort_values("padj", na_position="last")
    return out.drop_duplicates("gene").set_index("gene")
