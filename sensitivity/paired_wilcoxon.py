#!/usr/bin/env python
"""Paired sensitivity analysis: Wilcoxon signed-rank on tumour-normal pairs from the same patient.

The main analysis compares tumours and normals with the Mann-Whitney U (rank-sum) test, which
treats them as independent. In TCGA every solid-tissue normal comes from a patient who also has a
tumour in the data, so this script repeats the tests on those pairs only:

1. Design B differential methylation: all 38 COAD tumour-normal pairs, every probe tested in the
   main analysis (same matrix, same QC and imputation, same thresholds).
   Three calls per probe:
     main   Mann-Whitney U, 295 tumours vs 38 normals (the published call, from dmc_B.parquet)
     mwu38  Mann-Whitney U, the 38 paired tumours vs the 38 normals (same samples, pairing ignored)
     paired Wilcoxon signed-rank on the 38 within-patient differences
   mwu38 vs paired isolates the effect of pairing; main vs mwu38 the effect of dropping 257 tumours.
2. Ten-cancer specificity screen: the same tumour-normal comparison in each screened cancer,
   on its patient-matched pairs only.
3. The selection funnel rerun with the paired calls: promoter DMCs -> hyper-down genes ->
   candidate CpGs -> colon-specific CpGs -> blood-compatible CpGs -> 11-CpG panels.

Usage: paired_wilcoxon.py WORKDIR OUTDIR   (WORKDIR/cache holds the cache of a full run_all.py run)
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

HERE = Path(__file__).resolve().parent
# the pipeline sits beside this folder (crc_meth_repro/) or, in the public repository, one level up
PKG = HERE.parent / "crc_meth_repro" if (HERE.parent / "crc_meth_repro").is_dir() else HERE.parent
sys.path.insert(0, str(PKG / "src"))

from crcmeth import blood as B           # noqa: E402
from crcmeth import classify as C        # noqa: E402
from crcmeth import data as D            # noqa: E402
from crcmeth import integrate as I       # noqa: E402
from crcmeth import specificity as S     # noqa: E402
from crcmeth.dmc import bh               # noqa: E402

CHUNK = 20000


def load_cfg():
    import yaml
    return yaml.safe_load((PKG / "config.yaml").read_text())


def pairs(tumors, normals):
    """(tumour, normal) per patient with both; a patient with several tumours keeps the first."""
    t_by = {}
    for t in sorted(tumors):
        t_by.setdefault(D.patient_of(t), t)
    out = [(t_by[D.patient_of(n)], n) for n in sorted(normals) if D.patient_of(n) in t_by]
    pts = [D.patient_of(n) for _, n in out]
    assert len(pts) == len(set(pts)), "one normal per patient expected"
    return out


def signed_rank(diff: np.ndarray) -> np.ndarray:
    """Two-sided Wilcoxon signed-rank p per row. NaN pairs and zero differences are dropped
    (Wilcoxon's convention). Exact null distribution for up to 50 pairs without ties or zeros,
    otherwise the normal approximation with tie and continuity correction (scipy's own rule)."""
    p = np.full(diff.shape[0], np.nan)
    for i in range(0, diff.shape[0], CHUNK):
        d = diff[i:i + CHUNK]
        with np.errstate(all="ignore"):
            # scipy picks one method per call, so split rows: exact where |d| has no zeros or ties
            # and n <= 50, normal approximation otherwise; rows with NaN pairs handled one by one
            full = ~np.isnan(d).any(1)
            ad = np.sort(np.abs(d), axis=1)
            clean = full & (ad[:, 0] > 0) & (np.diff(ad, axis=1) > 0).all(1) & (d.shape[1] <= 50)
            r = np.full(len(d), np.nan)
            for rows, meth in ((clean, "exact"), (full & ~clean, "approx")):
                if rows.any():
                    r[rows] = wilcoxon(d[rows], axis=1, zero_method="wilcox", correction=True,
                                       method=meth).pvalue
            if (~full).any():
                r[~full] = wilcoxon(d[~full], axis=1, zero_method="wilcox", correction=True,
                                    method="auto", nan_policy="omit").pvalue
            p[i:i + CHUNK] = r
    return p


def rank_sum(t: np.ndarray, n: np.ndarray) -> np.ndarray:
    p = np.full(t.shape[0], np.nan)
    for i in range(0, t.shape[0], CHUNK):
        p[i:i + CHUNK] = mannwhitneyu(t[i:i + CHUNK], n[i:i + CHUNK], axis=1,
                                      method="asymptotic", nan_policy="omit").pvalue
    return p


def call(delta, fdr, cfg):
    d = cfg["dmc"]
    sig = (fdr < d["fdr"]) & (np.abs(delta) > d["min_abs_delta_beta"])
    return np.where(sig & (delta > 0), "hyper", np.where(sig, "hypo", "ns"))


# --------------------------------------------------------------------------- 1. Design B DMCs
def design_b(cache: Path, cfg: dict, coad: dict):
    beta = D.load_frame(cache / "coad_beta.parquet")            # QC'd + imputed discovery matrix
    main = D.load_frame(cache / "dmc_B.parquet")
    assert beta.index.equals(main.index), "probe set differs from the main DMC table"
    pr = pairs(coad["tumors"], coad["normals"])
    tcols, ncols = [t for t, _ in pr], [n for _, n in pr]
    t, n = beta[tcols].to_numpy(np.float64), beta[ncols].to_numpy(np.float64)
    diff = t - n
    out = pd.DataFrame(index=beta.index.rename("probe"))
    out["delta_main"] = main["delta"]
    out["p_main"] = main["p"]
    out["fdr_main"] = main["fdr"]
    out["status_main"] = main["status"]
    out["n_pairs"] = np.isfinite(diff).sum(1)
    out["delta_paired"] = np.nanmean(diff, 1)                # = mean(paired T) - mean(N)
    out["median_diff_paired"] = np.nanmedian(diff, 1)
    out["frac_pairs_tumor_higher"] = (diff > 0).sum(1) / out["n_pairs"]
    out["p_mwu38"] = rank_sum(t, n)
    out["fdr_mwu38"] = bh(out["p_mwu38"].to_numpy())
    out["status_mwu38"] = call(out["delta_paired"].to_numpy(), out["fdr_mwu38"].to_numpy(), cfg)
    out["p_wilcoxon"] = signed_rank(diff)
    out["fdr_wilcoxon"] = bh(out["p_wilcoxon"].to_numpy())
    out["status_paired"] = call(out["delta_paired"].to_numpy(), out["fdr_wilcoxon"].to_numpy(), cfg)
    # main-analysis columns the pipeline functions read, under the paired call
    dmc_paired = pd.DataFrame({"delta": out["delta_paired"], "fdr": out["fdr_wilcoxon"],
                               "p": out["p_wilcoxon"], "status": out["status_paired"],
                               "mean_t": np.nanmean(t, 1), "mean_n": np.nanmean(n, 1)},
                              index=beta.index)
    # the test alone swapped: paired FDR, Δβ of the main analysis (all 295 tumours)
    dmc_test = main.copy()
    dmc_test["fdr"] = out["fdr_wilcoxon"]
    dmc_test["p"] = out["p_wilcoxon"]
    dmc_test["status"] = call(main["delta"].to_numpy(), out["fdr_wilcoxon"].to_numpy(), cfg)
    out["status_paired_test_main_delta"] = dmc_test["status"]
    tumor_median = beta[coad["tumors"]].median(axis=1)      # blood rule uses all 295 tumours
    return out, dmc_paired, dmc_test, main, pr, tumor_median


# --------------------------------------------------------------------------- 2. screen
def paired_screen_stats(cache: Path, cfg: dict, probes) -> dict[str, pd.DataFrame]:
    """Per screened cancer, paired statistics for `probes` from the kept rows of the pan-cancer
    stream (raw betas, the same rows the screen reads). The raw matrices of the screened cancers
    were not kept, so paired p-values exist for the kept rows only (30,763 probes selected on the
    COAD tumour-normal difference). Two FDRs:
      q_kept  = BH over the kept rows. The rows were chosen on COAD data, independently of the
                screened cancer's own test, so BH over them controls the FDR (independent
                filtering); the family is smaller than the main screen's.
      q_upper = BH over the kept rows with m = every probe tested in that cancer by the main
                screen. It can only exceed the q-value of a genome-wide paired BH, so it gives
                the fewest exclusions (most conservative screen)."""
    pan = json.loads((cache / "pancan.done").read_text())["cancers"]
    max_na = cfg["probe_qc"]["max_missing_frac"]
    out = {}
    for c in cfg["screen_cancers"]:
        m = pan[c]
        kept = D.load_frame(cache / "pancan" / f"{c}.kept.parquet")
        st = D.load_frame(cache / "pancan" / f"{c}.stats.parquet")
        m_total = int(st["p"].notna().sum())
        pr = pairs(m["tumors"], m["normals"])
        t = kept[[a for a, _ in pr]].to_numpy(np.float64)
        n = kept[[b for _, b in pr]].to_numpy(np.float64)
        diff = t - n
        bad = np.isnan(diff).mean(1) > max_na
        p = np.full(len(kept), np.nan)
        if len(pr) >= 2:
            p[~bad] = signed_rank(diff[~bad])
        ok = np.isfinite(p)
        q_up, q_kept = np.full(len(kept), np.nan), np.full(len(kept), np.nan)
        if ok.any():
            q_kept[ok] = bh(p[ok])
            q_up[ok] = np.minimum(q_kept[ok] * m_total / ok.sum(), 1.0)   # BH scales linearly in m
        # the same paired samples, pairing ignored (rank-sum): separates pairing from sample choice
        p_mwu = np.full(len(kept), np.nan)
        if len(pr) >= 2:
            p_mwu[~bad] = rank_sum(t[~bad], n[~bad])
        q_mwu = np.full(len(kept), np.nan)
        okm = np.isfinite(p_mwu)
        if okm.any():
            q_mwu[okm] = bh(p_mwu[okm])
        with np.errstate(all="ignore"):
            delta = np.nanmean(diff, 1)
        df = pd.DataFrame({"n_pairs": np.isfinite(diff).sum(1), "delta_paired": delta,
                           "p_wilcoxon": p, "q_kept": q_kept, "q_upper": q_up, "p_mwu_pairs": p_mwu, "q_kept_mwu_pairs": q_mwu,
                           "delta_main_screen": st["delta"].reindex(kept.index),
                           "median_t": np.nanmedian(kept[m["tumors"]].to_numpy(np.float64), 1)},
                          index=kept.index)
        df = df.reindex(sorted(set(probes)))
        df.attrs.update(n_pairs=len(pr), m_total=m_total, m_kept=int(ok.sum()),
                        n_tumor=m["n_tumor"], n_normal=m["n_normal"],
                        min_p=float(min_two_sided_p(len(pr))))
        out[c] = df
    return out


def min_two_sided_p(n: int) -> float:
    """Smallest two-sided exact signed-rank p with n pairs: all signs equal, 2 / 2^n."""
    return min(1.0, 2.0 / 2 ** n) if n else 1.0


def paired_screen(cands, stats: dict, main_stats: dict, cfg: dict, qcol: str) -> pd.DataFrame:
    """S.screen with the paired call (q from column `qcol`). A cancer whose pairs cannot reach
    p < 0.05 even with every pair in the same direction (CESC 3 pairs, SARC 4) keeps its
    rank-sum call from the main screen, so the paired screen never loses a cancer by construction."""
    d = cfg["dmc"]
    out = pd.DataFrame(index=pd.Index(sorted(cands), name="probe"))
    for c, st in stats.items():
        st = st.reindex(out.index)
        if st.attrs.get("min_p", min_two_sided_p(stats[c].attrs["n_pairs"])) >= d["fdr"]:
            ms = main_stats[c].reindex(out.index)
            flag = ((ms["fdr"] < d["fdr"]) & (ms["delta"].abs() > d["min_abs_delta_beta"])).fillna(False)
            out[f"{c}_call"] = "rank-sum (too few pairs)"
        else:
            flag = ((st[qcol] < d["fdr"]) & (st["delta_paired"].abs() > d["min_abs_delta_beta"])).fillna(False)
            out[f"{c}_call"] = "signed-rank"
        out[f"{c}_excluded"] = flag
        out[f"{c}_delta_paired"] = st["delta_paired"]
        out[f"{c}_untestable"] = st["p_wilcoxon"].isna()
    flags = [c for c in out.columns if c.endswith("_excluded")]
    out["n_cancers_excluding"] = out[flags].sum(axis=1)
    out["specific"] = out["n_cancers_excluding"] == 0
    return out


# --------------------------------------------------------------------------- 3. funnel
def funnel(dmc, cache, cfg, deg_status, ann, screen_fn, bref, tumor_median):
    pdef = cfg["classify"]["panel_promoter_groups"]
    res = {}
    for name in ("primary", "extended"):
        groups = cfg["promoter_groups"][name]
        dmg = I.promoter_dmg(dmc, ann, groups)
        q = I.quadrants(dmg, deg_status)
        cands = I.candidate_cpgs(dmc, ann, q["hyper_down"], groups)
        scr = screen_fn(cands.index)
        spec = cands.loc[scr.index[scr["specific"]]]
        res[name] = {"n_dmc": int((dmc["status"] != "ns").sum()),
                     "n_hyper": int((dmc["status"] == "hyper").sum()),
                     "n_hypo": int((dmc["status"] == "hypo").sum()),
                     "hyper_genes": q["hyper_genes"], "hyper_down_genes": len(q["hyper_down"]),
                     "candidates": cands, "screen": scr, "specific": spec}
    spec = res[pdef]["specific"]
    have = spec.index.intersection(bref.index)
    margin = tumor_median.reindex(have) - bref.loc[have, "blood_beta"]
    ok = (bref.loc[have, "blood_beta"] <= cfg["blood"]["max_blood_beta"]) & (margin >= cfg["blood"]["min_margin"])
    size = cfg["classify"]["panel_size"]
    unknown = sorted(set(spec.index) - set(have))
    res["blood_known"] = len(have)
    res["blood_unknown"] = unknown
    res["blood_compatible"] = sorted(ok.index[ok])
    # blood data exist only for the probes of the main analysis; a CpG without them is bracketed:
    res["panel_blood_compatible"] = C.delta_panel(spec.loc[ok.index[ok]], size)          # unknown = not compatible
    res["panel_blood_if_unknown_ok"] = C.delta_panel(spec.loc[list(ok.index[ok]) + unknown], size)  # unknown = compatible
    res["unknown_in_panel_if_ok"] = sorted(set(res["panel_blood_if_unknown_ok"]) & set(unknown))
    res["panel_screened"] = C.delta_panel(spec, size)
    return res


def main(work, outdir):
    work, outdir = Path(work), Path(outdir)
    cache = work / "cache"
    outdir.mkdir(parents=True, exist_ok=True)
    cfg = load_cfg()
    coad = json.loads((cache / "coad.done").read_text())
    ann = pickle.loads((cache / "annot.pkl").read_bytes())
    deg = D.load_frame(cache / "deg.parquet")["status"]
    pdef = cfg["classify"]["panel_promoter_groups"]

    # ---- 1
    allcpg, dmc_paired, dmc_test, dmc_main, pr, tumor_median = design_b(cache, cfg, coad)
    print(f"Design B: {len(pr)} pairs, {len(allcpg)} probes")

    # blood reference and published panels
    # blood profile of every probe the main analysis fetched (B.blood_profile on the selection
    # half of GSE40279 and on GSE35069, as in run_all._blood_reference); checked against Table S6a
    wb, _, cells, cell_type, _, _ = pickle.loads((cache / "blood_inputs.pkl").read_bytes())
    order = D.rng(cfg, 41).permutation(wb.columns)
    bref = B.blood_profile(sorted(set(wb.index) | set(cells.index)), wb[order[: len(order) // 2]], cells, cell_type)
    s6 = pd.read_csv(PKG / "results" / "tables" / "blood_background_specific_B.csv", index_col=0)
    assert np.allclose(bref["blood_beta"].reindex(s6.index), s6["blood_beta"], equal_nan=True), "blood profile differs from S6a"
    bref = bref[bref["blood_beta"].notna()]
    panel_bc = json.loads((cache / "classify.done").read_text())["panels"]["blood_compatible_11"]
    panel_net = json.loads((cache / "network.done").read_text())["panel"]

    # ---- 2: paired screen stats for every CpG either funnel could screen
    cand_main_any = set()
    cand_paired_any = set()
    for name in ("primary", "extended"):
        g = cfg["promoter_groups"][name]
        cand_main_any |= set(I.candidate_cpgs(dmc_main, ann, I.quadrants(I.promoter_dmg(dmc_main, ann, g), deg)["hyper_down"], g).index)
        for dm in (dmc_paired, dmc_test):
            cand_paired_any |= set(I.candidate_cpgs(dm, ann, I.quadrants(I.promoter_dmg(dm, ann, g), deg)["hyper_down"], g).index)
    probes = cand_main_any | cand_paired_any | set(panel_bc) | set(panel_net)
    pstats = paired_screen_stats(cache, cfg, probes)
    pan_main = {c: D.load_frame(cache / "pancan" / f"{c}.stats.parquet") for c in cfg["screen_cancers"]}
    pan_meta = {c: m for c, m in json.loads((cache / "pancan.done").read_text())["cancers"].items()
                if c in cfg["screen_cancers"]}
    main_screen = lambda idx: S.screen(idx, pan_main, pan_meta, cfg)
    kept_screen = lambda idx: paired_screen(idx, pstats, pan_main, cfg, "q_kept")
    up_screen = lambda idx: paired_screen(idx, pstats, pan_main, cfg, "q_upper")
    mwu_screen = lambda idx: paired_screen(idx, pstats, pan_main, cfg, "q_kept_mwu_pairs")

    # ---- 3: funnels
    V = [("A", "main: rank-sum DMC (295 T vs 38 N), rank-sum screen", dmc_main, main_screen),
         ("B", "signed-rank FDR with the main Δβ (test swapped only), rank-sum screen", dmc_test, main_screen),
         ("C", "signed-rank FDR and paired Δβ (38 pairs), rank-sum screen", dmc_paired, main_screen),
         ("D", "as C, signed-rank screen (BH over the kept probes)", dmc_paired, kept_screen),
         ("E", "as C, signed-rank screen (genome-wide BH upper bound; fewest exclusions)", dmc_paired, up_screen),
         ("F", "as D, rank-sum on the same paired samples (pairing ignored)", dmc_paired, mwu_screen)]
    F = {f"{k}. {lab}": funnel(dm, cache, cfg, deg, ann, fn, bref, tumor_median) for k, lab, dm, fn in V}
    fm = F[next(iter(F))]
    assert fm["panel_blood_compatible"] == panel_bc, "main funnel does not reproduce the blood-compatible panel"
    assert len(fm[pdef]["specific"]) == 943, len(fm[pdef]["specific"])

    # ---------------------------------------------------------------- outputs
    allcpg.insert(0, "genes_extended", I.genes_in(ann.loc[allcpg.index], cfg["promoter_groups"]["extended"]).map(";".join))
    allcpg.to_csv(outdir / "S_paired_all_cpgs_designB.csv.gz", float_format="%.6g")

    # panels
    rows = []
    for pname, panel in (("blood-compatible", panel_bc), ("network-prioritized", panel_net)):
        for p in panel:
            r = allcpg.loc[p]
            rec = {"panel": pname, "probe": p, "genes": r["genes_extended"],
                   "delta_main_295v38": r["delta_main"], "p_main": r["p_main"], "fdr_main": r["fdr_main"], "status_main": r["status_main"],
                   "delta_paired_38": r["delta_paired"], "median_within_pair_diff": r["median_diff_paired"],
                   "pairs_tumor_higher": f"{int(round(r['frac_pairs_tumor_higher'] * r['n_pairs']))}/{int(r['n_pairs'])}",
                   "p_wilcoxon": r["p_wilcoxon"], "fdr_wilcoxon": r["fdr_wilcoxon"], "status_paired": r["status_paired"],
                   "fdr_mwu_same38": r["fdr_mwu38"],
                   "same_direction": np.sign(r["delta_main"]) == np.sign(r["delta_paired"]),
                   "dmc_in_both": r["status_main"] == r["status_paired"] == "hyper"}
            for c, st in pstats.items():
                rec[f"{c}_delta_paired"] = st.at[p, "delta_paired"]
                rec[f"{c}_q_kept"] = st.at[p, "q_kept"]
            for k in ("D", "E"):
                scr = [f for n_, f in F.items() if n_.startswith(k)][0][pdef]["screen"]
                rec[f"specific_{k}"] = bool(scr.at[p, "specific"]) if p in scr.index else "not a candidate"
                if p in scr.index:
                    rec[f"excluded_by_{k}"] = ";".join(c for c in pstats if scr.at[p, f"{c}_excluded"])
            for k, f in F.items():
                k = k[0]
                rec[f"candidate_{k}"] = p in f[pdef]["candidates"].index
                rec[f"in_blood_panel_{k}"] = p in f["panel_blood_compatible"]
            rows.append(rec)
    panels = pd.DataFrame(rows)

    # funnel summary
    summ = []
    for fname, f in F.items():
        e = f[pdef]
        spec = set(e["specific"].index)
        summ.append({"analysis": fname,
                     "DMCs": e["n_dmc"], "hyper DMCs": e["n_hyper"], "hypo DMCs": e["n_hypo"],
                     "hyper-down genes (core)": f["primary"]["hyper_down_genes"],
                     "hyper-down genes (extended)": e["hyper_down_genes"],
                     "candidate CpGs (extended)": len(e["candidates"]),
                     "candidates shared with main": len(set(e["candidates"].index) & set(fm[pdef]["candidates"].index)),
                     "colon-specific CpGs": len(spec),
                     "colon-specific shared with main": len(spec & set(fm[pdef]["specific"].index)),
                     **{f"excluded by {c}": int(e["screen"][f"{c}_excluded"].sum()) for c in cfg["screen_cancers"]},
                     "specific CpGs without blood data": len(f["blood_unknown"]),
                     "blood-compatible CpGs": len(f["blood_compatible"]),
                     "blood-compatible shared with main": len(set(f["blood_compatible"]) & set(fm["blood_compatible"])),
                     "blood-compatible panel": ";".join(f["panel_blood_compatible"]),
                     "panel overlap with main blood-compatible panel": len(set(f["panel_blood_compatible"]) & set(panel_bc)),
                     "panel if CpGs without blood data were compatible": ";".join(f["panel_blood_if_unknown_ok"]),
                     "its overlap with main panel": len(set(f["panel_blood_if_unknown_ok"]) & set(panel_bc)),
                     "screen-derived panel": ";".join(f["panel_screened"]),
                     "screen-derived panel overlap with main": len(set(f["panel_screened"]) & set(fm["panel_screened"])),
                     "network-prioritized CpGs still colon-specific": len(set(panel_net) & spec)})
    summ = pd.DataFrame(summ)

    # genome-wide concordance of the DMC calls
    st = allcpg
    conc = []
    for a, b in (("status_main", "status_paired"), ("status_mwu38", "status_paired"), ("status_main", "status_mwu38")):
        conc.append(pd.crosstab(st[a], st[b]).rename_axis(index=a, columns=b))
    dmc_main_set = st.index[st.status_main != "ns"]
    agree = {
        "probes": len(st),
        "pairs": len(pr),
        "DMCs main": int((st.status_main != "ns").sum()),
        "DMCs paired": int((st.status_paired != "ns").sum()),
        "DMCs mwu same 38 pairs (pairing ignored)": int((st.status_mwu38 != "ns").sum()),
        "DMCs signed-rank FDR with main delta": int((st.status_paired_test_main_delta != "ns").sum()),
        "main DMCs signed-rank FDR >= 0.05": int((st.loc[st.index[st.status_main != "ns"], "fdr_wilcoxon"] >= 0.05).sum()),
        "main DMCs kept by paired test": int((st.loc[dmc_main_set, "status_paired"] == st.loc[dmc_main_set, "status_main"]).sum()),
        "main DMCs with same sign of paired delta": int((np.sign(st.loc[dmc_main_set, "delta_paired"]) == np.sign(st.loc[dmc_main_set, "delta_main"])).sum()),
        "main DMCs paired FDR < 0.05": int((st.loc[dmc_main_set, "fdr_wilcoxon"] < 0.05).sum()),
        "main DMCs paired |delta| > 0.2": int((st.loc[dmc_main_set, "delta_paired"].abs() > 0.2).sum()),
        "Pearson r of delta (main vs paired), all probes": float(np.corrcoef(st.delta_main, st.delta_paired)[0, 1]),
        "screen: cancers and pairs": {c: s.attrs["n_pairs"] for c, s in pstats.items()},
    }

    # candidate CpGs that a paired screen cannot test (not in the kept rows of the stream)
    miss = {c: int(s.reindex(sorted(cand_paired_any))["p_wilcoxon"].isna().sum()) for c, s in pstats.items()}
    agree["paired candidates without a paired screen p (kept rows / missingness)"] = miss
    agree["screen: smallest attainable two-sided p"] = {c: s.attrs["min_p"] for c, s in pstats.items()}

    compact = panels[["panel", "probe", "genes", "delta_main_295v38", "p_main", "fdr_main",
                      "delta_paired_38", "median_within_pair_diff", "pairs_tumor_higher",
                      "p_wilcoxon", "fdr_wilcoxon", "status_paired", "same_direction"]].rename(columns={
        "probe": "CpG", "delta_main_295v38": "Δβ main (295 T − 38 N)", "p_main": "p main (rank-sum)",
        "fdr_main": "FDR main", "delta_paired_38": "Δβ paired (mean within-pair difference)",
        "median_within_pair_diff": "median within-pair difference", "pairs_tumor_higher": "pairs with T > N",
        "p_wilcoxon": "p signed-rank", "fdr_wilcoxon": "FDR signed-rank", "status_paired": "paired call",
        "same_direction": "same direction"})
    dm = allcpg[(allcpg.status_main != "ns") | (allcpg.status_paired != "ns")]
    # composition of the screened cancers: tumours with vs without a matched normal, main candidates
    comp = []
    cand_a = fm[pdef]["candidates"].index
    for c in cfg["screen_cancers"]:
        m = pan_meta[c]
        k = D.load_frame(cache / "pancan" / f"{c}.kept.parquet").reindex(cand_a)
        pt = {D.patient_of(x) for x in m["normals"]}
        tp = [x for x in m["tumors"] if D.patient_of(x) in pt]
        tu = [x for x in m["tumors"] if D.patient_of(x) not in pt]
        comp.append({"cancer": c, "tumours with a matched normal": len(tp), "tumours without": len(tu),
                     "median tumour β, matched tumours": float(np.nanmedian(k[tp].median(axis=1))) if tp else np.nan,
                     "median tumour β, unmatched tumours": float(np.nanmedian(k[tu].median(axis=1))) if tu else np.nan,
                     "median normal β": float(np.nanmedian(k[m["normals"]].median(axis=1))),
                     "candidates excluded, main screen": int(fm[pdef]["screen"][f"{c}_excluded"].sum()),
                     "candidates excluded, signed-rank on pairs (D)": int(F[[n_ for n_ in F if n_.startswith("D")][0]][pdef]["screen"][f"{c}_excluded"].sum()),
                     "candidates excluded, rank-sum on the same pairs (F)": int(F[[n_ for n_ in F if n_.startswith("F")][0]][pdef]["screen"][f"{c}_excluded"].sum())})
    comp = pd.DataFrame(comp)
    comp.to_csv(outdir / "screen_composition.csv", index=False)
    with pd.ExcelWriter(outdir / "paired_sensitivity.xlsx", engine="openpyxl") as xw:
        pd.DataFrame({"sheet": ["panels_22", "panels_detail", "funnel", "dmc_concordance", "screen_pairs",
                                "screen_composition", "panel_cpgs_screen", "all_dmcs"],
                      "content": [
                          "The 22 panel CpGs (11 blood-compatible, 11 network-prioritized; GRIN2A cg16368442 is in both): Δβ, p and FDR in the main analysis (Mann–Whitney U, 295 tumours vs 38 normals) and the paired analysis (Wilcoxon signed-rank, 38 tumour–normal pairs from the same patients).",
                          "Same CpGs with the rank-sum on the same 38 pairs (pairing ignored), per-cancer paired screen statistics and their status in each funnel variant.",
                          "Selection funnel (extended promoter definition) under variants A–F; counts, per-cancer exclusions, panels and their overlap with the main panels.",
                          "Genome-wide cross-tabulation of DMC calls (hyper / hypo / ns): main vs paired, rank-sum on the 38 pairs vs paired, main vs rank-sum on the 38 pairs.",
                          "Pairs per screened cancer, number of probes tested and the smallest p a signed-rank test can reach with that many pairs.",
                          "Screened cancers: methylation of the main candidate CpGs in tumours with and without a matched normal; exclusions per screen variant.",
                          "Per-cancer paired statistics of the 21 distinct panel CpGs.",
                          f"Every probe that is a DMC in the main or the paired analysis ({len(dm):,} rows); all {len(allcpg):,} probes are in S_paired_all_cpgs_designB.csv.gz."]}
                     ).to_excel(xw, sheet_name="README", index=False)
        compact.to_excel(xw, sheet_name="panels_22", index=False)
        panels.to_excel(xw, sheet_name="panels_detail", index=False)
        summ.to_excel(xw, sheet_name="funnel", index=False)
        r0 = 0
        for ct in conc:
            ct.to_excel(xw, sheet_name="dmc_concordance", startrow=r0)
            r0 += len(ct) + 4
        pd.DataFrame([{"cancer": c, "pairs": s.attrs["n_pairs"], "tumors": s.attrs["n_tumor"],
                       "normals": s.attrs["n_normal"], "probes tested (main screen)": s.attrs["m_total"],
                       "probes tested (paired screen, kept rows)": s.attrs["m_kept"],
                       "smallest attainable p (exact)": s.attrs["min_p"],
                       "paired screen call": "signed-rank" if s.attrs["min_p"] < 0.05 else "rank-sum kept (too few pairs)"}
                      for c, s in pstats.items()]).to_excel(xw, sheet_name="screen_pairs", index=False)
        pc = pd.concat({c: s.reindex(list(dict.fromkeys(panel_bc + panel_net))) for c, s in pstats.items()},
                       names=["cancer", "probe"])
        pc.to_excel(xw, sheet_name="panel_cpgs_screen")
        comp.to_excel(xw, sheet_name="screen_composition", index=False)
        dm.to_excel(xw, sheet_name="all_dmcs")
    panels.to_csv(outdir / "paired_panels.csv", index=False)
    pd.concat({c: s_ for c, s_ in pstats.items()}, names=["cancer", "probe"]).to_csv(
        outdir / "S_paired_screen_candidates.csv.gz", float_format="%.6g")
    for k, f in F.items():   # per-variant screen of the candidates (extended promoter definition)
        f[pdef]["screen"].to_csv(outdir / f"screen_variant_{k[0]}.csv.gz")
    summ.to_csv(outdir / "paired_funnel.csv", index=False)
    (outdir / "paired_summary.json").write_text(json.dumps(agree, indent=1, default=str))
    for fname, f in F.items():
        print(fname, "| cands", len(f[pdef]["candidates"]), "| specific", len(f[pdef]["specific"]),
              "| blood-compatible", len(f["blood_compatible"]), "| unknown blood", len(f["blood_unknown"]),
              "| overlap", len(set(f["panel_blood_compatible"]) & set(panel_bc)),
              len(set(f["panel_blood_if_unknown_ok"]) & set(panel_bc)), "| panel", f["panel_blood_compatible"],
              "| unknown in panel if ok", f["unknown_in_panel_if_ok"])
    print(json.dumps(agree, indent=1, default=str))
    print(panels[["panel", "probe", "genes", "delta_main_295v38", "delta_paired_38", "pairs_tumor_higher",
                  "fdr_main", "fdr_wilcoxon", "status_paired"]].to_string())


if __name__ == "__main__":
    main(*sys.argv[1:3])
