"""Panel classifiers, leakage-free cross-validation, external validation, benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedGroupKFold, train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------------------------------- models
class Oversampled(BaseEstimator, ClassifierMixin):
    """Random minority oversampling applied to the training rows only."""

    def __init__(self, estimator=None, random_state=0):
        self.estimator = estimator
        self.random_state = random_state

    def fit(self, X, y):
        y = np.asarray(y)
        rs = np.random.default_rng(self.random_state)
        classes, counts = np.unique(y, return_counts=True)
        idx = [np.flatnonzero(y == c) for c in classes]
        idx = np.concatenate([i if len(i) == counts.max() else
                              np.concatenate([i, rs.choice(i, counts.max() - len(i))])
                              for i in idx])
        self.estimator_ = clone(self.estimator).fit(np.asarray(X)[idx], y[idx])
        self.classes_ = self.estimator_.classes_
        return self

    def predict_proba(self, X):
        return self.estimator_.predict_proba(X)

    def predict(self, X):
        return self.estimator_.predict(X)


def make_model(name: str, cfg: dict, seed: int):
    """Impute (median, fitted in-fold) -> scale -> classifier."""
    c = cfg["classify"]
    if name == "logreg":
        clf = LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
    elif name == "mlp":   # ~ MATLAB patternnet(10): one hidden layer of 10 units
        clf = Oversampled(MLPClassifier(hidden_layer_sizes=tuple(c["mlp_hidden"]),
                                        max_iter=c["mlp_max_iter"], early_stopping=False,
                                        random_state=seed), random_state=seed)
    elif name == "rf":
        clf = RandomForestClassifier(n_estimators=c["rf_trees"], class_weight="balanced",
                                     n_jobs=-1, random_state=seed)
    else:
        raise ValueError(name)
    return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), clf)


# --------------------------------------------------------------------------- metrics
def wilson(k: int, n: int, z: float = norm.ppf(0.975)) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + z ** 2 / n
    mid = (p + z ** 2 / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return (float(mid - half), float(mid + half))


def bootstrap_auc(y, s, n_boot: int, rs: np.random.Generator) -> tuple[float, float]:
    y, s = np.asarray(y), np.asarray(s)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    if len(pos) == 0 or len(neg) == 0:
        return (np.nan, np.nan)
    aucs = [roc_auc_score(np.r_[np.ones(len(pos)), np.zeros(len(neg))],
                          np.r_[s[rs.choice(pos, len(pos))], s[rs.choice(neg, len(neg))]])
            for _ in range(n_boot)]           # stratified: keeps both classes in every resample
    return tuple(float(v) for v in np.percentile(aucs, [2.5, 97.5]))


def metrics(y, score, threshold: float = 0.5, n_boot: int = 0,
            rs: np.random.Generator | None = None) -> dict:
    """y: 1 = tumour (positive class), 0 = normal."""
    y = np.asarray(y).astype(int)
    score = np.asarray(score, dtype=float)
    pred = (score >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
    out = {"n": len(y), "n_pos": int(tp + fn), "n_neg": int(tn + fp),
           "tp": int(tp), "fn": int(fn), "tn": int(tn), "fp": int(fp),
           "sensitivity": tp / (tp + fn) if tp + fn else np.nan,
           "specificity": tn / (tn + fp) if tn + fp else np.nan,
           "accuracy": (tp + tn) / len(y),
           "auc": roc_auc_score(y, score) if 0 < y.sum() < len(y) else np.nan,
           "threshold": threshold}
    out["sensitivity_ci"] = wilson(tp, tp + fn)
    out["specificity_ci"] = wilson(tn, tn + fp)
    out["accuracy_ci"] = wilson(tp + tn, len(y))
    if n_boot:
        out["auc_ci"] = bootstrap_auc(y, score, n_boot, rs)
    return out


def youden_threshold(y, score) -> float:
    fpr, tpr, thr = roc_curve(y, score)
    return float(thr[np.argmax(tpr - fpr)])


# --------------------------------------------------------------------------- CV
# Every split is grouped by patient: a patient's tumour and normal are always on the same
# side, so no fold is tested on a patient it trained on. `groups` holds one patient ID per
# sample, aligned with y.
def group_folds(y, groups, n_splits: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Stratified k-fold over whole patients."""
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = list(cv.split(np.zeros(len(y)), y, groups))
    groups = np.asarray(groups)
    for tr, te in folds:
        assert not set(groups[tr]) & set(groups[te]), "a patient is on both sides of a split"
    return folds


def group_split(y, groups, train_frac: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Hold-out split of whole patients, stratified by each patient's set of labels
    (tumour only, normal only, or both), so both classes keep their proportions."""
    y, groups = np.asarray(y), np.asarray(groups)
    ids, inv = np.unique(groups, return_inverse=True)
    stratum = pd.Series(y).groupby(inv).agg(lambda v: ",".join(map(str, sorted(set(v)))))
    g_tr, _ = train_test_split(np.arange(len(ids)), train_size=train_frac,
                               stratify=stratum.to_numpy(), random_state=seed)
    is_tr = np.isin(inv, g_tr)
    return np.flatnonzero(is_tr), np.flatnonzero(~is_tr)


def locked_panel_cv(X: pd.DataFrame, y: np.ndarray, model: str, cfg: dict, seed: int, *, groups):
    """Repeated patient-grouped stratified k-fold with a fixed feature set; all
    preprocessing (imputation, scaling) in-fold.

    Returns per-repeat metrics and the out-of-fold scores of every repeat."""
    c = cfg["classify"]
    oof = np.full((c["cv_repeats"], len(y)), np.nan)
    for r in range(c["cv_repeats"]):
        for k, (tr, te) in enumerate(group_folds(y, groups, c["cv_folds"], seed + r)):
            m = make_model(model, cfg, seed + r * c["cv_folds"] + k).fit(X.iloc[tr], y[tr])
            oof[r, te] = m.predict_proba(X.iloc[te])[:, 1]
    rows = [metrics(y, oof[r]) for r in range(c["cv_repeats"])]
    return pd.DataFrame(rows), oof


def nested_cv(beta: pd.DataFrame, y: np.ndarray, select_panel, model: str, cfg: dict,
              seed: int, rs: np.random.Generator, *, groups):
    """Outer patient-grouped stratified folds; `select_panel(train_columns)` reruns the
    whole selection (QC and imputation, DMC, DEG, integration, specificity, panel rule)
    on the training samples only and returns a CpG list. beta: probes x samples, not
    imputed (the model pipeline imputes from the training fold)."""
    c = cfg["classify"]
    score = np.full(len(y), np.nan)
    panels = []
    cols = np.asarray(beta.columns)
    for k, (tr, te) in enumerate(group_folds(y, groups, c["nested_outer_folds"], seed)):
        assert not set(cols[tr]) & set(cols[te])
        panel = list(select_panel(list(cols[tr])))
        panels.append({"fold": k, "panel": panel, "n_train": len(tr), "n_test": len(te)})
        if not panel:
            continue
        Xtr, Xte = beta.loc[panel, cols[tr]].T, beta.loc[panel, cols[te]].T
        m = make_model(model, cfg, seed + k).fit(Xtr, y[tr])
        score[te] = m.predict_proba(Xte)[:, 1]
    ok = ~np.isnan(score)
    res = metrics(y[ok], score[ok], n_boot=c["n_boot_auc"], rs=rs)
    res["n_scored"] = int(ok.sum())
    return res, panels, score


def single_split(X: pd.DataFrame, y: np.ndarray, model: str, cfg: dict, seed: int, rs, *, groups):
    """Single patient-grouped stratified 70/30 hold-out, as a comparison with cross-validation."""
    tr, te = group_split(y, groups, cfg["classify"]["split_train_frac"], seed)
    m = make_model(model, cfg, seed).fit(X.iloc[tr], y[tr])
    return metrics(y[te], m.predict_proba(X.iloc[te])[:, 1],
                   n_boot=cfg["classify"]["n_boot_auc"], rs=rs)


# --------------------------------------------------------------------------- external
def frozen_external(X_train: pd.DataFrame, y_train, externals: dict, model: str,
                    cfg: dict, seed: int, rs) -> tuple[dict, object]:
    """Fit once on all TCGA, apply unchanged to every external set.

    externals: name -> (X, y); y may contain NaN for groups scored descriptively
    (e.g. adenomas), which get only the predicted-positive fraction."""
    m = make_model(model, cfg, seed).fit(X_train, y_train)
    out = {}
    for name, (X, y) in externals.items():
        s = m.predict_proba(X[X_train.columns])[:, 1]
        y = np.asarray(y, dtype=float)
        lab = ~np.isnan(y)
        res = metrics(y[lab], s[lab], n_boot=cfg["classify"]["n_boot_auc"], rs=rs) if lab.any() else {}
        if (~lab).any():
            k = int((s[~lab] >= 0.5).sum())
            res["unlabelled_n"] = int((~lab).sum())
            res["unlabelled_called_positive"] = k
            res["unlabelled_called_positive_ci"] = wilson(k, int((~lab).sum()))
        res["scores"] = s
        out[name] = res
    return out, m


# --------------------------------------------------------------------------- benchmarks
def quick_cv_auc(X: np.ndarray, y: np.ndarray, seed: int, folds: int = 5, *, groups) -> float:
    s = np.zeros(len(y))
    for tr, te in group_folds(y, groups, folds, seed):
        m = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(),
                          LogisticRegression(class_weight="balanced", max_iter=2000))
        s[te] = m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return roc_auc_score(y, s)


def null_panels(beta: pd.DataFrame, y: np.ndarray, pool, size: int, n: int,
                rs: np.random.Generator, seed: int, *, groups) -> np.ndarray:
    """CV AUC of `n` random panels of `size` CpGs drawn from `pool` (e.g. all hyper DMCs)."""
    pool = np.asarray(sorted(set(pool) & set(beta.index)))
    return np.array([quick_cv_auc(beta.loc[rs.choice(pool, size, replace=False)].T.to_numpy(),
                                  y, seed, groups=groups) for _ in range(n)])


def mean_beta_score(beta: pd.DataFrame, panel) -> pd.Series:
    """Training-free score: mean beta over the panel (samples as index)."""
    return beta.loc[list(panel)].mean(axis=0)


def delta_panel(candidates: pd.DataFrame, size: int = 11) -> list[str]:
    """Panel selection rule: the highest-delta CpG of each gene, genes in decreasing order
    of that delta, until `size` CpGs (ties broken by probe ID). A CpG annotated to several
    genes represents all of them.

    candidates: index probe, columns delta, genes (';'-joined)."""
    per_gene = (candidates.assign(gene=candidates["genes"].str.split(";")).explode("gene")
                .rename_axis("probe").reset_index()
                .sort_values(["delta", "probe"], ascending=[False, True]))
    panel, covered = [], set()
    for p, g in zip(per_gene["probe"], per_gene["gene"]):
        if g in covered or p in panel:
            continue
        panel.append(p)
        covered.update(candidates.at[p, "genes"].split(";"))   # a CpG counts for all its genes
        if len(panel) == size:
            break
    return panel


def network_panel(ranking: pd.DataFrame, candidates: pd.DataFrame, size: int = 11) -> list[str]:
    """Comparator rule: genes in decreasing PageRank order (ties by gene symbol), the highest-delta
    CpG of each gene (ties by probe ID), until `size` CpGs. A CpG annotated to several genes counts
    for all of them; a gene already covered by a chosen CpG is skipped.

    ranking: index gene, column pagerank. candidates: index probe, columns delta, genes (';'-joined)."""
    order = ranking.assign(g=ranking.index).sort_values(["pagerank", "g"], ascending=[False, True]).index
    by_gene = (candidates.assign(gene=candidates["genes"].str.split(";")).explode("gene")
               .rename_axis("probe").reset_index().sort_values(["delta", "probe"], ascending=[False, True]))
    panel, covered = [], set()
    for g in order:
        if g in covered:
            continue
        hit = by_gene.loc[by_gene["gene"] == g, "probe"]
        hit = [p for p in hit if p not in panel]
        if not hit:
            continue
        panel.append(hit[0])
        covered.update(candidates.at[hit[0], "genes"].split(";"))
        if len(panel) == size:
            break
    return panel
