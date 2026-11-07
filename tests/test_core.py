"""Unit tests on synthetic data; no network access needed."""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from crcmeth import blood as B  # noqa: E402
from crcmeth import classify as C  # noqa: E402
from crcmeth import data as D  # noqa: E402
from crcmeth import dmc as M  # noqa: E402
from crcmeth import integrate as I  # noqa: E402
from crcmeth import report as R  # noqa: E402
from crcmeth import specificity as S  # noqa: E402


@pytest.fixture(scope="module")
def cfg(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("run")
    return D.load_config(ROOT / "config.yaml", str(tmp / "w"), str(tmp / "r"))


# --------------------------------------------------------------------------- config / data
def test_config_round_trip(cfg):
    assert cfg["classify"]["panel_size"] == 11
    assert "network_prioritized_panel" not in cfg and cfg["network"]["alpha"] == 0.85
    assert len(cfg["screen_cancers"]) == 10
    assert cfg["designs"]["primary"] in ("A", "B")
    assert Path(cfg["workdir"]).is_dir() and Path(cfg["resultsdir"]).is_dir()


def test_barcodes():
    assert D.parse_barcode("TCGA-AA-3517-01A") == ("TCGA-AA-3517", 1, "A")
    cols = ["sample", "TCGA-AA-0001-01B", "TCGA-AA-0001-01A", "TCGA-AA-0001-11A",
            "TCGA-AA-0002-01A", "TCGA-AA-0003-02A", "garbage"]
    assert D.select_samples(cols, D.TUMOR_CODE) == ["TCGA-AA-0001-01A", "TCGA-AA-0002-01A"]
    assert D.select_samples(cols, D.NORMAL_CODE) == ["TCGA-AA-0001-11A"]


def test_series_matrix(tmp_path):
    txt = ("!Series_title\t\"x\"\n"
           "!Sample_geo_accession\t\"GSM1\"\t\"GSM2\"\n"
           "!Sample_characteristics_ch1\t\"disease status: cancer\"\t\"disease status: normal-H\"\n"
           "!series_matrix_table_begin\n"
           "\"ID_REF\"\t\"GSM1\"\t\"GSM2\"\n"
           "cg1\t0.1\t0.2\ncg2\t0.9\tnull\n"
           "!series_matrix_table_end\n")
    p = tmp_path / "s.txt.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(txt)
    md, tab = D.parse_series_matrix(p)
    assert md.loc["GSM2", "disease status"] == "normal-H"
    assert tab.shape == (2, 2) and np.isnan(tab.loc["cg2", "GSM2"])
    sub = D.series_table(p, ["cg2", "cg9"], chunk_rows=1)
    assert list(sub.index) == ["cg2"] and sub.dtypes.iloc[0] == np.float32


# --------------------------------------------------------------------------- statistics
def test_bh_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = np.random.default_rng(0).random(500) ** 3
    assert np.allclose(M.bh(p), multipletests(p, method="fdr_bh")[1])
    p[[3, 7]] = np.nan
    q = M.bh(p)
    assert np.isnan(q[[3, 7]]).all()
    ok = ~np.isnan(p)
    assert np.allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1])


def test_dmc_table_matches_scipy(cfg):
    rs = np.random.default_rng(1)
    t = rs.beta(2, 5, (300, 40))
    n = rs.beta(2, 5, (300, 12))
    t[:30] += 0.4                                  # 30 true hyper probes
    beta = pd.DataFrame(np.c_[t, n].clip(0, 1), index=[f"cg{i}" for i in range(300)],
                        columns=[f"T{i}" for i in range(40)] + [f"N{i}" for i in range(12)])
    tumor, normal = list(beta.columns[:40]), list(beta.columns[40:])
    dt = M.dmc_table(beta, tumor, normal, {**cfg, "dmc": {**cfg["dmc"], "chunk_rows": 64}})
    for i in (0, 77, 299):
        ref = mannwhitneyu(beta.iloc[i][tumor], beta.iloc[i][normal], method="asymptotic").pvalue
        assert dt["p"].iloc[i] == pytest.approx(ref)
    assert (dt["status"].iloc[:30] == "hyper").all()
    assert (dt["status"].iloc[30:] != "hyper").mean() > 0.95
    assert np.allclose(dt["delta"], beta[tumor].mean(1) - beta[normal].mean(1))


def test_impute_is_label_free():
    b = pd.DataFrame([[0.1, np.nan, 0.3], [np.nan, 0.5, 0.7]], columns=list("abc"))
    out = M.impute_rows(b)
    assert out.loc[0, "b"] == pytest.approx(0.2) and out.loc[1, "a"] == pytest.approx(0.6)
    # permuting columns (labels) does not change the filled values
    assert np.allclose(M.impute_rows(b[["c", "a", "b"]])[list("abc")], out)


def test_probe_qc(cfg):
    ann = pd.DataFrame({"chrom": ["1", "X", "2", "3"]}, index=["cg1", "cg2", "cg3", "cg4"])
    b = pd.DataFrame({"s1": [0.1, 0.2, 0.3, np.nan], "s2": [0.1, 0.2, 0.3, 0.4]},
                     index=["cg1", "cg2", "cg3", "cg4"])
    out = M.probe_qc(b, ann, {"cg3"}, cfg)
    assert list(out.index) == ["cg1"]              # X dropped, cg3 masked, cg4 50% missing
    # sample-independent filters only: cg4 is left for the training fold to judge
    assert list(M.probe_qc(b, ann, {"cg3"}, cfg, missing=False).index) == ["cg1", "cg4"]


def test_filter_impute_fits_on_given_columns(cfg):
    b = pd.DataFrame({"tr1": [0.1, np.nan, np.nan], "tr2": [0.3, 0.5, np.nan],
                      "tr3": [np.nan, 0.7, 0.2], "te1": [0.9, 0.9, 0.9]}, index=["a", "b", "c"])
    train = ["tr1", "tr2", "tr3"]
    c = {**cfg, "probe_qc": {**cfg["probe_qc"], "max_missing_frac": 0.5}}
    out = M.filter_impute(b[train], c)
    assert list(out.index) == ["a", "b"]           # c: 2/3 missing in training
    assert out.loc["a", "tr3"] == pytest.approx(0.2) and out.loc["b", "tr1"] == pytest.approx(0.6)
    # fitted on all columns, the test column would change both the probes kept and the fills
    leaky = M.filter_impute(b, c)
    assert "c" in leaky.index and leaky.loc["a", "tr3"] == pytest.approx(0.3)


# --------------------------------------------------------------------------- integration
def _toy_annot():
    return pd.DataFrame({
        "island": ["Island", "N_Shore", "OpenSea", "Island", "Island"],
        "gene_group": [[("G1", "TSS1500")], [("G1", "TSS200")], [("G2", "Body")],
                       [("G3", "TSS200")], [("G4", "1stExon"), ("MIR1", "TSS1500")]]},
        index=["c1", "c2", "c3", "c4", "c5"])


def test_promoter_integration():
    ann = _toy_annot()
    dmc = pd.DataFrame({"status": ["hyper", "hyper", "hyper", "hypo", "hyper"],
                        "delta": [0.3, 0.4, 0.5, -0.3, 0.25], "fdr": 1e-5,
                        "mean_t": 0.6, "mean_n": 0.2}, index=ann.index)
    dmg = I.promoter_dmg(dmc, ann, ["TSS1500", "TSS200"])
    assert dmg.loc["G1", "n_hyper"] == 2 and dmg.loc["G3", "status"] == "hypo"
    assert "G2" not in dmg.index and "G4" not in dmg.index
    deg = pd.Series({"G1": "down", "G3": "up", "MIR1": "down"})
    q = I.quadrants(dmg, deg)
    assert q["hyper_down"] == ["G1", "MIR1"] and q["hypo_up"] == ["G3"]
    cand = I.candidate_cpgs(dmc, ann, q["hyper_down"], ["TSS1500", "TSS200"])
    assert list(cand.index) == ["c1", "c2", "c5"] and cand.loc["c5", "genes"] == "MIR1"
    ext = I.candidate_cpgs(dmc, ann, ["G4"], ["TSS1500", "TSS200", "1stExon"])
    assert list(ext.index) == ["c5"]


def test_specificity_screen(cfg):
    st = pd.DataFrame({"delta": [0.3, 0.05, 0.25], "fdr": [1e-4, 1e-4, 0.2],
                       "median_t": [0.6, 0.5, 0.1]}, index=["c1", "c2", "c3"])
    meta = {"X": {"n_normal": 10}, "Y": {"n_normal": 3}}
    out = S.screen(["c1", "c2", "c3"], {"X": st, "Y": st * 0 + 0.0}, meta, cfg)
    assert out["specific"].to_dict() == {"c1": False, "c2": True, "c3": True}
    assert out.attrs["underpowered"] == ["Y"]
    strict = S.screen(["c1", "c2", "c3"], {"X": st}, meta,
                      {**cfg, "specificity": {**cfg["specificity"], "rule": "dmc_or_high"}})
    assert strict["specific"].to_dict() == {"c1": False, "c2": False, "c3": True}


# --------------------------------------------------------------------------- classification
def _toy_xy(n_t=60, n_n=20, p=5, seed=0):
    rs = np.random.default_rng(seed)
    y = np.r_[np.ones(n_t), np.zeros(n_n)].astype(int)
    X = rs.normal(size=(n_t + n_n, p)) + y[:, None] * 1.5
    return pd.DataFrame(X, columns=[f"cg{i}" for i in range(p)],
                        index=[f"s{i}" for i in range(n_t + n_n)]), y


def test_metrics_and_wilson():
    from statsmodels.stats.proportion import proportion_confint
    y = np.array([1, 1, 1, 0, 0, 0, 0, 1])
    s = np.array([0.9, 0.8, 0.2, 0.1, 0.6, 0.3, 0.2, 0.7])
    m = C.metrics(y, s)
    assert (m["tp"], m["fn"], m["tn"], m["fp"]) == (3, 1, 3, 1)
    assert m["sensitivity"] == m["specificity"] == m["accuracy"] == 0.75
    assert np.allclose(C.wilson(3, 4), proportion_confint(3, 4, method="wilson"))


def test_oversampled_balances():
    from sklearn.dummy import DummyClassifier
    X, y = _toy_xy()
    m = C.Oversampled(DummyClassifier(strategy="prior"), random_state=0).fit(X.values, y)
    assert np.allclose(m.predict_proba(X.values[:1]), [[0.5, 0.5]])


def _toy_patients(n_t=60, n_n=20):
    """As in TCGA: tumour i and normal i come from patient i, so the normals are all paired."""
    return np.array([f"P{i}" for i in range(n_t)] + [f"P{i}" for i in range(n_n)])


def test_locked_cv_deterministic_and_sane(cfg):
    X, y = _toy_xy()
    g = _toy_patients()
    c = {**cfg, "classify": {**cfg["classify"], "cv_repeats": 2, "cv_folds": 5}}
    df1, oof1 = C.locked_panel_cv(X, y, "logreg", c, 3, groups=g)
    df2, oof2 = C.locked_panel_cv(X, y, "logreg", c, 3, groups=g)
    assert np.array_equal(oof1, oof2)
    assert not np.isnan(oof1).any() and df1["auc"].min() > 0.8
    assert not np.array_equal(oof1[0], oof1[1])        # repeats use different folds


def test_splits_keep_each_patient_on_one_side():
    _, y = _toy_xy()
    g = _toy_patients()
    for seed in range(5):
        folds = C.group_folds(y, g, 5, seed)
        assert sorted(np.concatenate([te for _, te in folds]).tolist()) == list(range(len(y)))
        for tr, te in folds:
            assert not set(g[tr]) & set(g[te])
            assert 0 < (y[te] == 0).sum() and 0 < (y[te] == 1).sum()
        tr, te = C.group_split(y, g, 0.7, seed)
        assert not set(g[tr]) & set(g[te]) and len(tr) + len(te) == len(y)
        assert abs((y[tr] == 0).mean() - (y == 0).mean()) < 0.05      # stratified
    # sample-level stratified k-fold would split pairs; the grouped one never does
    from sklearn.model_selection import StratifiedKFold
    split = [set(g[tr]) & set(g[te]) for tr, te in
             StratifiedKFold(5, shuffle=True, random_state=0).split(y, y)]
    assert any(split)


def test_nested_cv_never_sees_test_samples(cfg):
    X, y = _toy_xy()
    g = _toy_patients()
    beta = X.T
    seen = []

    def select(train_cols):
        seen.append(set(train_cols))
        return ["cg0", "cg1"]

    c = {**cfg, "classify": {**cfg["classify"], "nested_outer_folds": 4, "n_boot_auc": 50}}
    res, panels, score = C.nested_cv(beta, y, select, "logreg", c, 0, np.random.default_rng(0), groups=g)
    assert len(seen) == 4
    patient = dict(zip(beta.columns, g))
    for k, tr in enumerate(seen):
        test = set(beta.columns) - tr
        assert len(test) == panels[k]["n_test"] and not (tr & test)
        assert not {patient[s] for s in tr} & {patient[s] for s in test}
    # every sample is scored exactly once, by a model that did not train on it
    assert res["n_scored"] == len(y) and not np.isnan(score).any()


def test_model_pipeline_imputes_from_training_rows():
    X, y = _toy_xy()
    g = _toy_patients()
    X = X.copy()
    X.iloc[::7, 0] = np.nan                       # missing values reach the model unimputed
    assert C.quick_cv_auc(X.to_numpy(), y, 0, groups=g) > 0.8
    tr = np.r_[0:30, 60:70]                       # 30 tumours and 10 normals
    m = C.make_model("logreg", {"classify": {}}, 0).fit(X.iloc[tr], y[tr])
    fill = m.named_steps["simpleimputer"].statistics_[0]
    assert fill == pytest.approx(X.iloc[tr, 0].median()) != X.iloc[:, 0].median()


def test_frozen_external_unlabelled_group(cfg):
    X, y = _toy_xy()
    Xe, ye = _toy_xy(10, 10, seed=5)
    ye = ye.astype(float)
    ye[:3] = np.nan                              # e.g. adenomas
    c = {**cfg, "classify": {**cfg["classify"], "n_boot_auc": 50}}
    res, _ = C.frozen_external(X, y, {"ext": (Xe, ye)}, "logreg", c, 0, np.random.default_rng(0))
    assert res["ext"]["n"] == 17 and res["ext"]["unlabelled_n"] == 3


def test_delta_panel_one_cpg_per_gene():
    cand = pd.DataFrame({"delta": [0.3, 0.5, 0.4, 0.2, 0.4], "genes": ["A", "A", "B;C", "C", "D"]},
                        index=["p1", "p2", "p3", "p4", "p5"])
    assert C.delta_panel(cand, 3) == ["p2", "p3", "p5"]      # p3 covers B and C; p3 < p5 on the tie
    assert C.delta_panel(cand, 10) == ["p2", "p3", "p5"]     # fewer genes than the requested size


# --------------------------------------------------------------------------- blood background
def test_blood_profile_and_compatibility(cfg):
    wb = pd.DataFrame({"s1": [0.02, 0.60], "s2": [0.04, 0.70], "s3": [0.05, 0.65]}, index=["lo", "hi"])
    cells = pd.DataFrame({"c1": [0.03, 0.8], "c2": [0.30, 0.9]}, index=["lo", "hi"])
    prof = B.blood_profile(["lo", "hi", "absent"], wb, cells, pd.Series({"c1": "T", "c2": "Mono"}))
    assert prof.loc["lo", "celltype_argmax"] == "Mono"
    # one leukocyte type methylating the CpG is enough to disqualify it
    assert prof.loc["lo", "blood_beta"] == pytest.approx(0.30)
    ok = B.blood_compatible(prof, pd.Series({"lo": 0.6, "hi": 0.95, "absent": 0.9}), cfg)
    assert not ok.any()
    prof.loc["lo", "blood_beta"] = 0.05
    ok = B.blood_compatible(prof, pd.Series({"lo": 0.6, "hi": 0.95, "absent": 0.9}), cfg)
    assert ok.tolist() == [True, False, False]


def test_dilution_monotone_and_specificity_fixed():
    rs = np.random.default_rng(0)
    idx = [f"p{i}" for i in range(5)]
    tumor = pd.DataFrame(rs.uniform(0.6, 0.9, (5, 40)), index=idx)
    blood = pd.DataFrame(rs.uniform(0.0, 0.05, (5, 200)), index=idx)
    sim = B.simulate_dilution(tumor, blood.iloc[:, :100], blood.iloc[:, 100:], {"P": idx},
                              [0.001, 0.05, 0.2], [0, 100], 400, 0.95, np.random.default_rng(1))
    for _, g in sim.groupby("depth"):
        assert g["sensitivity_at_spec"].is_monotonic_increasing
        assert g["auc"].iloc[-1] > 0.99
    lod = B.detection_limit(sim, 0.8)
    assert set(lod["lowest_fraction_detected"]) <= {0.05, 0.2}


def test_moss_groups_and_cfdna_scores():
    md = pd.DataFrame({
        "title": ["cfDNA_healthy_individual_1", "cfDNA_young_male_pool_rep1", "cfDNA_colon_cancer_806",
                  "cfDNA_colon_cancer_806_50ng", "cfDNA_mix_5", "in_vitro_mix_12", "in_vitro_mix_9",
                  "cfDNA_sepsis_2", "colon_epithelial_cells_rep1", "leukocytes_1"],
        "source_name_ch1": ["cfDNA, normal", "cfDNA of young males, pooled from multiple donors",
                            "cfDNA, colon adenocarcinoma", "cfDNA, colon adenocarcinoma, 50 ng DNA used",
                            "cfDNA, in vitro mix", "in vitro mix of genomic DNA", "in vitro mix of genomic DNA",
                            "cfDNA, sepsis", "colon epithelial cells", "leukocytes"],
        "description": ["HI1", "Young males", "Metastases: Liver, lung, bone", "Metastases: Liver, lung, bone",
                        "10% pl806, 90% HI2", "10% Colon cells, 5% Neurons", "10% Hepatocytes, 3.5% Lung cells",
                        "AST:29", "Colon 1", "HI2_leuk"]},
        index=[f"GSM{i}" for i in range(10)])
    g = B.moss_groups(md)
    assert g["group"].tolist() == ["healthy_cfdna", "healthy_cfdna", "colon_cancer_cfdna", "colon_cancer_cfdna",
                                   "cfdna_spike_in", "gdna_mix", "gdna_mix", "sepsis_cfdna",
                                   "colon_epithelium", "leukocytes"]
    assert g.loc["GSM3", "replicate_of"] == "cfDNA_colon_cancer_806" and g.loc["GSM2", "replicate_of"] == ""
    assert g.loc["GSM4", "spike_fraction"] == pytest.approx(0.10)
    assert g.loc["GSM5", "colon_cell_fraction"] == pytest.approx(0.10)
    assert g.loc["GSM6", "colon_cell_fraction"] == 0.0
    beta = pd.DataFrame(np.full((2, 10), 0.02), index=["a", "b"], columns=g.index)
    beta.loc[:, "GSM1"] = 0.04
    beta.loc[:, "GSM2"] = 0.50
    sc = B.cfdna_scores(beta, g, {"P": ["a", "b", "not_on_array"]}).set_index("sample")
    assert sc["n_cpgs_on_array"].iloc[0] == 2
    assert sc.loc["GSM2", "above_all_healthy"] and not sc.loc["GSM7", "above_all_healthy"]


class _FigCtx:
    """Minimal stand-in for run_all.Ctx: cache and results under one temporary directory."""

    def __init__(self, cfg, root):
        import pickle
        self._pickle = pickle
        self.cfg, self.res, self.cache = cfg, root / "res", root / "cache"
        for p in (self.res / "figures", self.res / "tables", self.cache):
            p.mkdir(parents=True, exist_ok=True)

    def path(self, name):
        return self.cache / name

    def dump(self, name, obj):
        self.path(name).write_bytes(self._pickle.dumps(obj))

    def load(self, name):
        return self._pickle.loads(self.path(name).read_bytes())


def test_new_figures_render(cfg, tmp_path, monkeypatch):
    try:
        import pyarrow  # noqa: F401
    except ImportError:                       # no parquet engine here: frames go through pickle
        monkeypatch.setattr(D, "save_frame", lambda df, path: df.to_pickle(path))
        monkeypatch.setattr(D, "load_frame", pd.read_pickle)
    ctx = _FigCtx(cfg, tmp_path)
    d, pdef = cfg["designs"]["primary"], cfg["classify"]["panel_promoter_groups"]
    ann = _toy_annot()
    ctx.dump("annot.pkl", ann)
    dmc = pd.DataFrame({"status": ["hyper", "hyper", "hyper", "hypo", "hyper"],
                        "delta": [0.3, 0.4, 0.5, -0.3, 0.25], "fdr": 1e-5}, index=ann.index)
    D.save_frame(dmc, ctx.path(f"dmc_{d}.parquet"))
    D.save_frame(pd.DataFrame({"log2fc": [-2.0, 2.0, -1.8], "status": ["down", "up", "down"]},
                              index=["G1", "G3", "MIR1"]), ctx.path("deg.parquet"))
    dmg = I.promoter_dmg(dmc, ann, cfg["promoter_groups"][pdef])
    ctx.dump(f"integrate_{d}_{pdef}.pkl",
             {"quadrants": I.quadrants(dmg, pd.Series({"G1": "down", "G3": "up", "MIR1": "down"}))})
    pd.DataFrame({"Term": ["Axonogenesis (GO:0007409)", "Wnt signaling pathway"],
                  "Overlap": ["5/100", "3/80"], "Adjusted P-value": [1e-4, 0.01],
                  "library": ["GO_Biological_Process_2023", "KEGG_2021_Human"]}
                 ).to_csv(ctx.res / "tables" / "enrichment_hyper_down.csv", index=False)
    R.fig_integration(ctx)
    R.fig_enrichment(ctx)
    for name in ("fig_integration", "fig_enrichment"):
        assert (ctx.res / "figures" / f"{name}.png").stat().st_size > 10_000
    assert R._name("network_prioritized_11") == "network-prioritized panel" and R._name("other") == "other"


def test_network_panel_rule():
    """Genes in PageRank order, best-delta CpG each; a CpG shared by two genes covers both."""
    cand = pd.DataFrame({"delta": [0.30, 0.25, 0.22, 0.21, 0.40],
                         "genes": ["A", "A", "B;C", "D", "E"]},
                        index=["p1", "p2", "p3", "p4", "p5"])
    rk = pd.DataFrame({"pagerank": [0.5, 0.4, 0.3, 0.2, 0.1]}, index=["A", "C", "B", "D", "E"])
    assert C.network_panel(rk, cand, 3) == ["p1", "p3", "p4"]      # C covers B; E never reached
    assert C.network_panel(rk, cand, 11) == ["p1", "p3", "p4", "p5"]
