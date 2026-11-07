"""Traceability of the comparator panel, figures and the run summary."""

from __future__ import annotations


import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_curve  # noqa: E402

from . import data as D  # noqa: E402
from .annot import genes_in  # noqa: E402

def traceability(ctx, panel) -> pd.DataFrame:
    """Status of each comparator CpG at every pipeline stage."""
    cfg = ctx.cfg
    ann, mask = ctx.load("annot.pkl"), ctx.load("mask.pkl")
    rows = []
    for p in panel:
        r = {"probe": p, "annotation": ";".join(f"{g}:{x}" for g, x in ann.at[p, "gene_group"]),
             "island": ann.at[p, "island"], "zhou_masked": p in mask}
        for d in ("A", "B"):
            dmc = D.load_frame(ctx.path(f"dmc_{d}.parquet"))
            r[f"{d}_passes_qc"] = p in dmc.index
            if p in dmc.index:
                r[f"{d}_delta"] = round(float(dmc.at[p, "delta"]), 4)
                r[f"{d}_fdr"] = float(dmc.at[p, "fdr"])
                r[f"{d}_status"] = dmc.at[p, "status"]
            for pdef in ("primary", "extended"):
                cand = ctx.load(f"integrate_{d}_{pdef}.pkl")["candidates"]
                r[f"{d}_{pdef}_candidate"] = p in cand.index
                scr = ctx.load(f"screen_{d}_{pdef}.pkl")
                if p in scr.index:
                    r[f"{d}_{pdef}_specific"] = bool(scr.at[p, "specific"])
                    r[f"{d}_{pdef}_excluded_by"] = ";".join(
                        c.replace("_excluded", "") for c in scr.columns
                        if c.endswith("_excluded") and scr.at[p, c])
        rows.append(r)
    return pd.DataFrame(rows).set_index("probe")


# --------------------------------------------------------------------------- figures
# Panel names as they appear in the manuscript
PANEL_NAMES = {
    "blood_compatible_11": "blood-compatible panel",
    "screened_11": "screen-derived panel",
    "network_prioritized_11": "network-prioritized panel",
    "network_prioritized_blood_compatible": "blood-compatible CpGs of the network-prioritized panel (3)",
    "SEPT9_promoter": "SEPT9 promoter",
}
PANEL_COLORS = {"network_prioritized_11": "#8e44ad", "screened_11": "#2ca02c", "blood_compatible_11": "#d62728",
                "SEPT9_promoter": "#000000", "network_prioritized_blood_compatible": "#1f77b4"}
MODEL_NAMES = {"logreg": "logistic regression", "mlp": "neural network", "rf": "random forest"}


def _name(panel: str) -> str:
    return PANEL_NAMES.get(panel, panel)


def _save(fig, ctx, name, dpi=300):
    fig.tight_layout()
    fig.savefig(ctx.res / "figures" / f"{name}.png", dpi=dpi)
    plt.close(fig)


REGION_LABELS = {"CGI": "CpG island", "Island": "island", "CGI_promoter": "CGI promoter", "N_Shelf": "N shelf",
                 "N_Shore": "N shore", "S_Shelf": "S shelf", "S_Shore": "S shore", "OpenSea": "open sea",
                 "1stExon": "1st exon", "Body": "gene body"}


def fig_regions(ctx, info):
    reg = info["dmc"][ctx.cfg["designs"]["primary"]]["regions"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, key, title in zip(axes, ["fig1A_all", "fig1B_island", "fig1C_region"],
                              ["A  DMC class", "B  CpG-island relation", "C  Gene region"]):
        if key == "fig1A_all":
            d = {k.replace("fig1A_", ""): reg[k] for k in reg if k.startswith("fig1A_")}
        else:
            d = reg[key]
        df = pd.DataFrame(d).T.fillna(0).reindex(columns=["hyper", "hypo"], fill_value=0)
        df.index = [REGION_LABELS.get(i, i) for i in df.index]
        df.plot.bar(ax=ax, color=["#c0392b", "#2c7fb8"], width=0.8)
        ax.set_title(title)
        ax.set_ylabel("CpGs")
        ax.tick_params(axis="x", rotation=45)
    _save(fig, ctx, "fig1_dmc_regions")


def fig_volcano(ctx):
    d = ctx.cfg["designs"]["primary"]
    dmc = D.load_frame(ctx.path(f"dmc_{d}.parquet"))
    fig, ax = plt.subplots(figsize=(5, 4))
    y = -np.log10(dmc["fdr"].clip(lower=1e-300))
    ax.hexbin(dmc["delta"], y, gridsize=120, bins="log", cmap="Greys", mincnt=1)
    pub = dmc.index.intersection(ctx.load_info("network")["panel"])
    ax.scatter(dmc.loc[pub, "delta"], y[pub], c="#c0392b", s=18, label="network-prioritized panel (11 CpGs)")
    ax.axvline(0.2, ls="--", c="k", lw=0.7)
    ax.axvline(-0.2, ls="--", c="k", lw=0.7)
    ax.set_xlabel("Δβ (tumor − normal)")
    ax.set_ylabel("−log10 FDR")
    ax.set_title(f"TCGA-COAD, design {d}")
    ax.legend(loc="upper left", fontsize=8)
    _save(fig, ctx, "fig2_volcano")


QUAD_COLORS = {"hyper_down": "#c0392b", "hyper_up": "#e6a6a0", "hypo_up": "#2c7fb8", "hypo_down": "#a6c8e0"}


def fig_integration(ctx):
    """A: Δβ of promoter DMCs against log2 fold change of their gene; B: genes per quadrant."""
    cfg = ctx.cfg
    d, pdef = cfg["designs"]["primary"], cfg["classify"]["panel_promoter_groups"]
    dmc = D.load_frame(ctx.path(f"dmc_{d}.parquet"))
    deg = D.load_frame(ctx.path("deg.parquet"))
    q = ctx.load(f"integrate_{d}_{pdef}.pkl")["quadrants"]
    dm = dmc[dmc["status"] != "ns"]
    pg = genes_in(ctx.load("annot.pkl").loc[dm.index], cfg["promoter_groups"][pdef]).explode().dropna()
    pairs = pd.DataFrame({"probe": pg.index, "gene": pg.values})
    pairs = pairs[pairs["gene"].isin(deg.index)]
    pairs["delta"] = dm["delta"].reindex(pairs["probe"]).to_numpy()
    pairs["log2fc"] = deg["log2fc"].reindex(pairs["gene"]).to_numpy()
    quad = {g: k for k in QUAD_COLORS for g in q[k]}
    pairs["quadrant"] = pairs["gene"].map(quad)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    rest = pairs[pairs["quadrant"].isna()]
    ax.scatter(rest["delta"], rest["log2fc"], s=3, color="0.8", rasterized=True, label="other pairs")
    for k, c in QUAD_COLORS.items():
        sub = pairs[pairs["quadrant"] == k]
        ax.scatter(sub["delta"], sub["log2fc"], s=5, color=c, rasterized=True,
                   label=f"{k.replace('_', '–')} ({len(sub):,} pairs)")
    lfc = cfg["deg"]["min_abs_log2fc"]
    for v in (-lfc, lfc):
        ax.axhline(v, ls="--", c="k", lw=0.6)
    ax.axvline(0, c="k", lw=0.6)
    ax.set_xlabel("Δβ (tumor − normal), promoter DMC")
    ax.set_ylabel("log2 fold change of the gene")
    ax.set_title("A  Promoter CpG–gene pairs")
    ax.legend(fontsize=7, markerscale=3, loc="upper right")
    ax = axes[1]
    x = np.arange(2)
    up = [len(q["hyper_up"]), len(q["hypo_up"])]
    down = [len(q["hyper_down"]), len(q["hypo_down"])]
    b1 = ax.bar(x - 0.2, up, 0.4, color=["#e6a6a0", "#2c7fb8"], edgecolor="k", label="upregulated")
    b2 = ax.bar(x + 0.2, down, 0.4, color=["#c0392b", "#a6c8e0"], edgecolor="k", hatch="//",
                label="downregulated")
    ax.bar_label(b1, fontsize=8)
    ax.bar_label(b2, fontsize=8)
    ax.set_xticks(x, [f"hypermethylated\n({q['hyper_genes']:,} genes)",
                      f"hypomethylated\n({q['hypo_genes']:,} genes)"])
    ax.set_ylabel("differentially expressed genes")
    ax.set_title("B  Genes per methylation–expression group")
    ax.legend(fontsize=8)
    _save(fig, ctx, "fig_integration")


def fig_enrichment(ctx, n=10):
    """Leading GO biological processes and KEGG pathways of the hyper-down genes."""
    enr = pd.read_csv(ctx.res / "tables" / "enrichment_hyper_down.csv")
    libs = [("GO_Biological_Process_2023", "GO biological process"), ("KEGG_2021_Human", "KEGG pathway")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, (lib, title) in zip(axes, libs):
        e = enr[enr["library"] == lib].nsmallest(n, "Adjusted P-value").iloc[::-1]
        y = -np.log10(e["Adjusted P-value"])
        k = e["Overlap"].str.split("/").str[0].astype(int)
        ax.barh(range(len(e)), y, color="#c0392b" if lib.startswith("GO") else "#2c3e50")
        ax.set_yticks(range(len(e)), [t.split(" (GO:")[0] for t in e["Term"]], fontsize=8)
        for i, (v, c) in enumerate(zip(y, k)):
            ax.text(v, i, f" {c}", va="center", fontsize=7)
        ax.axvline(-np.log10(0.05), ls="--", c="k", lw=0.6)
        ax.set_xlabel("−log10 FDR (number of genes at bar end)")
        ax.set_title(title)
    _save(fig, ctx, "fig_enrichment")


def fig_pancancer(ctx):
    cfg = ctx.cfg
    panels = ctx.load_info("classify")["panels"]
    dmc = D.load_frame(ctx.path(f"dmc_{cfg['designs']['primary']}.parquet"))
    stats = {"COAD": dmc["delta"]}
    for c in cfg["screen_cancers"]:
        stats[c] = D.load_frame(ctx.path("pancan") / f"{c}.stats.parquet")["delta"]
    keys = [p for p in ("blood_compatible_11", "network_prioritized_11") if p in panels]
    fig, axes = plt.subplots(1, len(keys), figsize=(6.5 * len(keys), 4.5), squeeze=False)
    for ax, p in zip(axes[0], keys):
        df = pd.DataFrame({c: st.reindex(panels[p]) for c, st in stats.items()}).T
        im = ax.imshow(df.values, cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
        ax.set_xticks(range(df.shape[1]), df.columns, rotation=60, fontsize=7)
        ax.set_yticks(range(df.shape[0]), df.index)
        ax.set_title(_name(p)[0].upper() + _name(p)[1:])
    fig.colorbar(im, ax=list(axes[0]), label="Δβ, tumor − own normal")
    fig.savefig(ctx.res / "figures" / "fig4_pancancer_delta.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_roc(ctx):
    scores = ctx.load("roc_scores.pkl")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, ds in zip(axes, ["TCGA_oof", "GSE48684", "GSE131013"]):
        for (panel, model, d), (y, s) in scores.items():
            if d != ds or model not in ("logreg", "mlp"):
                continue
            if len(np.unique(y)) < 2:
                continue
            fpr, tpr, _ = roc_curve(y, s)
            ax.plot(fpr, tpr, lw=1.2, color=PANEL_COLORS.get(panel), ls="-" if model == "logreg" else "--",
                    label=f"{_name(panel)}, {MODEL_NAMES.get(model, model)}")
        ax.plot([0, 1], [0, 1], "k:", lw=0.7)
        ax.set_title(ds.replace("_oof", " (out-of-fold)"))
        ax.set_xlabel("1 − specificity")
        ax.set_ylabel("sensitivity")
        ax.legend(fontsize=6, loc="lower right")
    _save(fig, ctx, "fig6_roc")


def fig_null(ctx, info):
    null = pd.read_csv(ctx.res / "tables" / "null_panel_auc.csv")["null_auc"]
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.hist(null, bins=40, color="#bdc3c7")
    for p, o in info["classify"]["null"]["observed"].items():
        ax.axvline(o["auc"], color=PANEL_COLORS.get(p), lw=1.5, label=f"{_name(p)}: {o['auc']:.3f} (p={o['empirical_p']:.3g})")
    ax.set_xlabel("5-fold CV AUC (logistic regression)")
    ax.set_ylabel("random 11-CpG hyper-DMC panels")
    ax.legend(fontsize=7)
    _save(fig, ctx, "fig7_null_panels")


CFDNA_LABELS = {"healthy_cfdna": "healthy", "colon_cancer_cfdna": "colon cancer",
                "cfdna_spike_in": "cfDNA spike-in", "lung_cancer_cfdna": "lung cancer",
                "breast_cancer_cfdna": "breast cancer", "cup_cfdna": "CUP", "sepsis_cfdna": "sepsis"}


def fig_blood(ctx):
    t = ctx.res / "tables"
    bg = pd.read_csv(t / "blood_background_panels.csv", index_col=0)
    sim = pd.read_csv(t / "blood_dilution.csv")
    cf = pd.read_csv(t / "cfdna_scores.csv")
    if "replicate_of" in cf.columns:  # repeated hybridizations are excluded from all group summaries
        cf = cf[cf["replicate_of"].isna()]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    pub = pd.concat([bg[bg[f"in_{p}"]] for p in ("blood_compatible_11", "network_prioritized_11")
                     if f"in_{p}" in bg.columns])
    x = np.arange(len(pub))
    ax.bar(x - 0.2, pub["tcga_tumor_median"], 0.4, label="TCGA tumor median", color="#c0392b")
    ax.bar(x + 0.2, pub["blood_beta"], 0.4, label="β$_{blood}$", color="#2c3e50")
    ax.axhline(ctx.cfg["blood"]["max_blood_beta"], color="k", ls=":", lw=0.8)
    if "in_blood_compatible_11" in bg.columns:
        ax.axvline(int(bg["in_blood_compatible_11"].sum()) - 0.5, color="grey", lw=0.8)
    ax.set_xticks(x, [f"{p}\n{g}" for p, g in zip(pub.index, pub["genes"].fillna(""))], rotation=90, fontsize=6)
    ax.set_ylabel("β")
    ax.set_title("A  Tumor and blood β (left: blood-compatible; right: network-prioritized)", fontsize=9)
    ax.legend(fontsize=7)
    ax = axes[1]
    for p, g in sim[sim["depth"] == 0].groupby("panel"):
        ax.plot(g["tumor_fraction"], g["sensitivity_at_spec"], "o-", color=PANEL_COLORS.get(p), label=_name(p))
    for p, g in sim[sim["depth"] > 0].groupby("panel"):
        ax.plot(g["tumor_fraction"], g["sensitivity_at_spec"], "o--", color=PANEL_COLORS.get(p), alpha=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("tumor DNA fraction in blood (solid: array; dashed: 100 reads)")
    ax.set_ylabel(f"sensitivity at {100 * ctx.cfg['blood']['specificity']:.0f}% specificity")
    ax.set_title("B  In-silico dilution into whole blood")
    ax.legend(fontsize=6)
    ax = axes[2]
    order = ["healthy_cfdna", "colon_cancer_cfdna", "cfdna_spike_in", "lung_cancer_cfdna",
             "breast_cancer_cfdna", "cup_cfdna", "sepsis_cfdna"]
    panels = [p for p in ("network_prioritized_11", "blood_compatible_11") if p in set(cf["panel"])]
    for j, p in enumerate(panels):
        d = cf[cf["panel"] == p]
        for i, g in enumerate(order):
            v = d.loc[d["group"] == g, "score"]
            ax.scatter(np.full(len(v), i + 0.2 * (j - 0.5)), v, s=14, color=PANEL_COLORS.get(p),
                       label=_name(p) if i == 0 else None)
    ax.set_xticks(range(len(order)), [CFDNA_LABELS.get(o, o) for o in order],
                  rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("mean β over panel")
    ax.set_title("C  Plasma cfDNA (GSE122126)")
    ax.legend(fontsize=7)
    _save(fig, ctx, "fig8_blood_cfdna")


# --------------------------------------------------------------------------- build
def build(ctx) -> dict:
    stages = ["annot", "coad", "pancan", "dmc", "rna", "integrate", "specificity", "network", "geo", "classify"]
    if ctx.path("blood.done").exists():
        stages.append("blood")
    info = {s: ctx.load_info(s) for s in stages}
    cv = pd.read_csv(ctx.res / "tables" / "classification_tcga.csv")
    ext = pd.read_csv(ctx.res / "tables" / "classification_external.csv")
    tr = traceability(ctx, ctx.load_info("network")["panel"])
    ctx.table(tr, "traceability_network_prioritized_panel.csv")
    errors = {}
    figs = {"regions": lambda: fig_regions(ctx, info), "volcano": lambda: fig_volcano(ctx),
            "integration": lambda: fig_integration(ctx), "enrichment": lambda: fig_enrichment(ctx),
            "pancancer": lambda: fig_pancancer(ctx),
            "roc": lambda: fig_roc(ctx), "null": lambda: fig_null(ctx, info)}
    if "blood" in info:
        figs["blood"] = lambda: fig_blood(ctx)
    for name, f in figs.items():
        try:
            f()
        except Exception as exc:          # one broken figure must not hide the tables
            errors[name] = repr(exc)
    D.write_json(info, ctx.res / "metrics.json")
    summary = _summary_md(ctx, info, cv, ext, tr)
    (ctx.res / "SUMMARY.md").write_text(summary)
    return {"figure_errors": errors, "n_tables": len(list((ctx.res / "tables").glob("*.csv")))}


def _fmt(x, pct=False):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "–"
    if pct:
        return f"{100 * x:.1f}"
    return f"{x:.3f}" if isinstance(x, float) else str(x)


def _summary_md(ctx, info, cv, ext, tr) -> str:
    L = ["# Analysis summary", ""]
    L += ["## TCGA performance (tumour = positive class)", ""]
    cols = ["panel", "model", "eval", "n", "sensitivity", "specificity", "accuracy", "auc", "auc_ci"]
    L += [cv[[c for c in cols if c in cv.columns]].to_markdown(index=False, floatfmt=".3f"), ""]
    L += ["## External validation, frozen TCGA model", ""]
    cols = ["panel", "model", "dataset", "subset", "n_pos", "n_neg", "tp", "fn", "tn", "fp",
            "sensitivity", "specificity", "accuracy", "auc", "unlabelled_called_positive", "unlabelled_n"]
    L += [ext[[c for c in cols if c in ext.columns]].to_markdown(index=False, floatfmt=".3f"), ""]
    n = info["classify"]["nested"]
    L += ["## Nested CV (whole selection funnel redone inside each outer fold)", ""]
    for label, block, key in (("with blood step", n, "network_prioritized_overlap"),
                              ("without blood step", n["screened"], "network_prioritized_overlap_screened")):
        for m in ("logreg", "mlp"):
            r = block[m]
            L.append(f"- {label}, {m}: AUC {_fmt(r['auc'])} (95% CI {', '.join(_fmt(v) for v in r['auc_ci'])}), "
                     f"sens {_fmt(r['sensitivity'], True)}%, spec {_fmt(r['specificity'], True)}%, "
                     f"comparator-CpG overlap per fold {[f[key] for f in r['folds']]}")
    nu = info["classify"]["null"]
    L += ["", "## Random-panel benchmark", "",
          f"Null median AUC {nu['null_auc_median']:.3f} (5–95%: {nu['null_auc_q05_q95'][0]:.3f}–"
          f"{nu['null_auc_q05_q95'][1]:.3f}) over {ctx.cfg['classify']['n_null_panels']} random panels "
          f"from {nu['n_pool']} hyper DMCs."]
    for p, o in nu["observed"].items():
        L.append(f"- {p}: AUC {o['auc']:.3f}, empirical p = {o['empirical_p']:.3g}")
    L += ["## Traceability of the 11 comparator CpGs", "", tr.to_markdown(), ""]
    L += [f"Determinism check (two identical runs): {info['classify']['deterministic']}", ""]
    if "blood" in info:
        L += _blood_md(ctx, info["blood"])
    return "\n".join(L)


def _blood_md(ctx, bl) -> list:
    t = ctx.res / "tables"
    b = ctx.cfg["blood"]
    L = ["## Blood background and plasma cfDNA", "",
         f"Whole blood GSE40279 (n = {bl['n_whole_blood']}; {bl['n_whole_blood_selection']} set the rule, "
         f"the rest are simulated), sorted leukocytes GSE35069 (n = {bl['n_sorted_cells']}), "
         f"plasma cfDNA GSE122126 EPIC (n = {bl['n_epic_samples']}). Blood-compatible: blood β ≤ "
         f"{b['max_blood_beta']} and tumour median − blood β ≥ {b['min_margin']}.", "",
         f"- CRC-specific CpGs passing: {bl['n_specific_blood_compatible']}/{bl['n_specific']}",
         f"- comparator CpGs passing: {len(bl['network_prioritized_blood_compatible'])}/11 "
         f"({', '.join(bl['network_prioritized_blood_compatible']) or 'none'})",
         f"- comparator CpGs absent from EPIC: {', '.join(bl['network_prioritized_not_on_epic']) or 'none'}", ""]
    bg = pd.read_csv(t / "blood_background_panels.csv", index_col=0)
    cols = ["genes", "tcga_tumor_median", "tcga_normal_median", "wb_median", "wb_p95",
            "celltype_max_median", "celltype_argmax", "blood_beta", "margin", "blood_compatible",
            "moss_healthy_cfdna_median", "moss_colon_cancer_cfdna_median"]
    L += ["### Network-prioritized 11 CpGs", "",
          bg.loc[bg["in_network_prioritized_11"], [c for c in cols if c in bg.columns]].to_markdown(floatfmt=".3f"), ""]
    lod = pd.read_csv(t / "blood_detection_limit.csv")
    L += [f"### In-silico dilution: lowest tumour fraction with ≥{100 * b['lod_sensitivity']:.0f}% sensitivity "
          f"at {100 * b['specificity']:.0f}% specificity", "", lod.to_markdown(index=False), ""]
    sim = pd.read_csv(t / "blood_dilution.csv")
    L += ["### Dilution grid", "", sim.drop(columns=["specificity"]).to_markdown(index=False, floatfmt=".3f"), ""]
    cf = pd.read_csv(t / "cfdna_scores.csv")
    keep = cf["group"].isin(["healthy_cfdna", "colon_cancer_cfdna", "cfdna_spike_in", "gdna_mix"]) | cf["group"].str.endswith("cancer_cfdna")
    L += ["### Real plasma cfDNA (z vs healthy cfDNA)", "",
          cf.loc[keep, ["panel", "n_cpgs_on_array", "title", "group", "spike_fraction", "colon_cell_fraction",
                        "score", "z_vs_healthy", "above_all_healthy"]].to_markdown(index=False, floatfmt=".3f"), ""]
    grp = (cf.groupby(["panel", "group"])["score"].agg(["count", "median", "min", "max"]).reset_index())
    L += ["### cfDNA score by group", "", grp.to_markdown(index=False, floatfmt=".3f"), ""]
    cls = pd.read_csv(t / "blood_panel_classification.csv")
    cols = ["panel", "model", "dataset", "n", "sensitivity", "specificity", "accuracy", "auc"]
    L += ["### Blood-compatible panels in tissue", "",
          cls[[c for c in cols if c in cls.columns]].to_markdown(index=False, floatfmt=".3f"), "",
          f"Panels: {bl['panels'].get('blood_compatible_11')}", ""]
    return L

