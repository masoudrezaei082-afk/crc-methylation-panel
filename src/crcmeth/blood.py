"""Blood background of tissue markers: leukocyte methylation, in-silico cfDNA dilution, real cfDNA.

Plasma cell-free DNA is mostly leukocyte-derived, so a tissue marker is usable in blood only
when leukocytes leave it unmethylated. This module scores each CpG against whole blood and
sorted leukocytes, simulates tumour DNA diluted into blood, and scores real plasma cfDNA."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


# --------------------------------------------------------------------------- background
def blood_profile(probes, whole_blood: pd.DataFrame, cells: pd.DataFrame,
                  cell_type: pd.Series) -> pd.DataFrame:
    """Per-CpG blood methylation. whole_blood / cells: probes x samples; cell_type: sample -> type."""
    probes = list(probes)
    wb = whole_blood.reindex(probes)
    ct = cells.reindex(probes).T.groupby(cell_type.reindex(cells.columns).values).median().T
    out = pd.DataFrame({"wb_median": wb.median(axis=1), "wb_p95": wb.quantile(0.95, axis=1),
                        "wb_n": wb.notna().sum(axis=1)}, index=probes)
    out["celltype_max_median"] = ct.max(axis=1)
    out["celltype_argmax"] = ct.fillna(-1.0).idxmax(axis=1).where(ct.notna().any(axis=1))
    out["blood_beta"] = out[["wb_p95", "celltype_max_median"]].max(axis=1)
    return out.join(ct.add_prefix("median_"))


def blood_compatible(profile: pd.DataFrame, tumor_median: pd.Series, cfg: dict) -> pd.Series:
    """A CpG suits blood when leukocytes leave it unmethylated and tumours methylate it well above that."""
    b = cfg["blood"]
    margin = tumor_median.reindex(profile.index) - profile["blood_beta"]
    return ((profile["blood_beta"] <= b["max_blood_beta"]) & (margin >= b["min_margin"])).rename("blood_compatible")


# --------------------------------------------------------------------------- dilution
def _observe(beta: np.ndarray, depth: int, rs: np.random.Generator) -> np.ndarray:
    """Array-like readout (depth 0) or bisulfite read counts at `depth` reads per CpG."""
    beta = np.clip(beta, 0.0, 1.0)
    if depth <= 0:
        return beta
    return rs.binomial(depth, beta) / depth


def simulate_dilution(tumor: pd.DataFrame, blood_cases: pd.DataFrame, blood_controls: pd.DataFrame,
                      panels: dict, fractions, depths, n: int, spec: float,
                      rs: np.random.Generator) -> pd.DataFrame:
    """Score = mean beta over the panel. Cases mix a random tumour into a random blood sample at
    tumour fraction f; controls are other individuals' blood. Sensitivity is read at the control
    quantile `spec`, so specificity is fixed by construction.

    All frames are probes x samples; missing values are filled with the probe median."""
    def fill(df):
        return df.T.fillna(df.median(axis=1)).T
    tumor, blood_cases, blood_controls = fill(tumor), fill(blood_cases), fill(blood_controls)
    rows = []
    for pname, panel in panels.items():
        panel = [p for p in panel if p in tumor.index and p in blood_cases.index
                 and p in blood_controls.index]
        if not panel:
            continue
        T = tumor.loc[panel].to_numpy()
        B = blood_cases.loc[panel].to_numpy()
        K = blood_controls.loc[panel].to_numpy()
        ti = rs.integers(0, T.shape[1], n)
        bi = rs.integers(0, B.shape[1], n)
        ki = rs.integers(0, K.shape[1], n)
        for depth in depths:
            ctrl = _observe(K[:, ki], depth, rs).mean(axis=0)
            thr = float(np.quantile(ctrl, spec))
            for f in fractions:
                case = _observe(f * T[:, ti] + (1 - f) * B[:, bi], depth, rs).mean(axis=0)
                y = np.r_[np.ones(n), np.zeros(n)]
                rows.append({"panel": pname, "n_cpgs": len(panel), "depth": depth, "tumor_fraction": f,
                             "auc": float(roc_auc_score(y, np.r_[case, ctrl])),
                             "sensitivity_at_spec": float((case > thr).mean()),
                             "specificity": spec, "threshold": thr,
                             "case_mean": float(case.mean()), "control_mean": float(ctrl.mean())})
    return pd.DataFrame(rows)


def detection_limit(sim: pd.DataFrame, min_sens: float) -> pd.DataFrame:
    """Lowest simulated tumour fraction reaching `min_sens` sensitivity, per panel and depth."""
    ok = sim[sim["sensitivity_at_spec"] >= min_sens]
    lod = ok.groupby(["panel", "depth"])["tumor_fraction"].min()
    return (sim[["panel", "depth"]].drop_duplicates().set_index(["panel", "depth"])
            .join(lod.rename("lowest_fraction_detected")).reset_index())


# --------------------------------------------------------------------------- real cfDNA
def moss_groups(md: pd.DataFrame) -> pd.DataFrame:
    """Sample groups of GSE122126 (Moss et al. 2018) from titles, source names and descriptions."""
    title = md["title"].str.lower()
    src = md["source_name_ch1"].str.lower()
    desc = md.get("description", pd.Series("", index=md.index)).fillna("")
    rules = [("healthy_cfdna", title.str.startswith("cfdna_healthy") | title.str.contains("_pool_")),
             ("colon_cancer_cfdna", src.str.startswith("cfdna, colon adenocarcinoma")),
             ("lung_cancer_cfdna", src.str.startswith("cfdna, lung")),
             ("breast_cancer_cfdna", src.str.startswith("cfdna, breast")),
             ("cup_cfdna", src.str.startswith("cfdna, cup")),
             ("sepsis_cfdna", src.str.contains("sepsis")),
             ("islet_transplant_cfdna", src.str.contains("islet transplant")),
             ("cfdna_spike_in", src.str.contains("cfdna, in vitro mix")),
             ("gdna_mix", src.str.startswith("in vitro mix")),
             ("colon_epithelium", src.str.startswith("colon epithelial")),
             ("leukocytes", src.str.startswith("leukocytes"))]
    group = pd.Series("other", index=md.index)
    for name, hit in reversed(rules):
        group[hit.values] = name
    spike = desc.map(lambda d: _pct(d, r"(\d+(?:\.\d+)?)% pl806"))
    colon = desc.map(lambda d: _pct(d, r"(\d+(?:\.\d+)?)% colon cells"))
    out = pd.DataFrame({"title": md["title"], "group": group, "description": desc,
                        "spike_fraction": spike, "colon_cell_fraction": colon})
    out.loc[out["group"] == "gdna_mix", "colon_cell_fraction"] = out["colon_cell_fraction"].fillna(0.0)
    # lower-input or repeated hybridisations of one plasma sample are not independent donors
    rep = out["title"].str.contains(r"_(?:50|100)ng|ng_rep2$", regex=True) & out["group"].str.endswith("cfdna")
    out["replicate_of"] = np.where(rep, out["title"].str.replace(r"_(?:50|100)ng.*|_rep2$", "", regex=True), "")
    return out


def _pct(text: str, pattern: str) -> float:
    m = re.search(pattern, text, flags=re.I)
    return float(m.group(1)) / 100 if m else np.nan


def cfdna_scores(beta: pd.DataFrame, groups: pd.DataFrame, panels: dict,
                 reference: str = "healthy_cfdna") -> pd.DataFrame:
    """Mean-beta panel score per sample and its z-score against the reference group."""
    rows = []
    ref = groups.index[groups["group"] == reference]
    for pname, panel in panels.items():
        panel = [p for p in panel if p in beta.index]
        if not panel:
            continue
        s = beta.loc[panel].mean(axis=0)
        mu, sd, top = s[ref].mean(), s[ref].std(ddof=1), s[ref].max()
        for smp, v in s.items():
            rows.append({"panel": pname, "n_cpgs_on_array": len(panel), "sample": smp,
                         **groups.loc[smp, ["title", "group", "spike_fraction",
                                            "colon_cell_fraction", "replicate_of"]].to_dict(),
                         "score": float(v), "z_vs_healthy": float((v - mu) / sd) if sd > 0 else np.nan,
                         "above_all_healthy": bool(v > top)})
    return pd.DataFrame(rows)
