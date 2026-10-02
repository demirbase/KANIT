"""Metrics from out-of-fold predictions.

An out-of-fold table has one row per genome per repeat of one arm, with columns
``repeat``, ``fold``, ``genome_id``, ``y`` (1 = resistant) and ``p`` (predicted
probability of resistance). The primary metric is the pooled out-of-fold ROC-AUC:
the AUC over all predictions of a repeat, averaged over the repeats. Confidence
intervals come from lineage-cluster bootstrap resamples in which the primary
metric is computed again; the same resamples serve both arms, so the interval of
their difference is paired.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)


def _auc(y, p, w=None) -> float:
    y = np.asarray(y)
    if w is not None:
        present = np.asarray(w) > 0
        if len(set(y[present].tolist())) < 2:
            return float("nan")
        return float(roc_auc_score(y[present], np.asarray(p)[present],
                                   sample_weight=np.asarray(w)[present]))
    if len(set(y.tolist())) < 2:
        return float("nan")
    return float(roc_auc_score(y, p))


def per_repeat(oof: pd.DataFrame, threshold: float = 0.5) -> pd.DataFrame:
    """Pooled metrics of every repeat."""
    rows = []
    for r, d in oof.groupby("repeat", sort=True):
        y, p = d["y"].to_numpy(), d["p"].to_numpy()
        pred = (p >= threshold).astype(int)
        rows.append({"repeat": r, "n": len(d), "roc_auc": _auc(y, p),
                     "pr_auc": float(average_precision_score(y, p)),
                     "pr_auc_baseline": float(y.mean()),
                     "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
                     "mcc": float(matthews_corrcoef(y, pred)),
                     "brier": float(brier_score_loss(y, p))})
    return pd.DataFrame(rows)


def per_fold(oof: pd.DataFrame) -> pd.DataFrame:
    """Fold-level AUC (undefined for a single-class fold) and fold sizes."""
    rows = []
    for (r, k), d in oof.groupby(["repeat", "fold"], sort=True):
        rows.append({"repeat": r, "fold": k, "n": len(d), "n_resistant": int(d["y"].sum()),
                     "roc_auc": _auc(d["y"], d["p"])})
    return pd.DataFrame(rows)


def summary(oof: pd.DataFrame, threshold: float = 0.5) -> dict:
    rep, fold = per_repeat(oof, threshold), per_fold(oof)
    return {
        "roc_auc": float(rep["roc_auc"].mean()),
        "roc_auc_sd_across_repeats": float(rep["roc_auc"].std(ddof=1)),
        "fold_roc_auc_mean": float(fold["roc_auc"].mean()),
        "fold_roc_auc_sd": float(fold["roc_auc"].std(ddof=1)),
        "n_folds_single_class": int(fold["roc_auc"].isna().sum()),
        "pr_auc": float(rep["pr_auc"].mean()),
        "pr_auc_baseline": float(rep["pr_auc_baseline"].mean()),
        "balanced_accuracy": float(rep["balanced_accuracy"].mean()),
        "mcc": float(rep["mcc"].mean()),
        "brier": float(rep["brier"].mean()),
        "threshold": threshold,
        "n_repeats": int(len(rep)),
    }


def reliability(oof: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    """Reliability curve of every repeat: equal-width bins of the predicted probability."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    rows = []
    for r, d in oof.groupby("repeat", sort=True):
        b = np.clip(np.digitize(d["p"], edges[1:-1]), 0, n_bins - 1)
        for k in range(n_bins):
            m = b == k
            rows.append({"repeat": r, "bin": k, "bin_low": edges[k], "bin_high": edges[k + 1],
                         "n": int(m.sum()),
                         "mean_predicted": float(d["p"][m].mean()) if m.any() else float("nan"),
                         "observed_rate": float(d["y"][m].mean()) if m.any() else float("nan")})
    return pd.DataFrame(rows)


def cluster_bootstrap(oofs: dict[str, pd.DataFrame], lineage_of: dict, *,
                      n_boot: int = 1000, seed: int = 0, diff: tuple[str, str] | None = None,
                      level: float = 0.95) -> dict:
    """Percentile intervals of the primary metric of every arm, and of
    ``diff[0] - diff[1]``, from lineage-cluster bootstrap resamples.

    Lineages are drawn with replacement; a genome's weight is the number of
    times its lineage was drawn. Every arm sees the same draws.
    """
    lineages = sorted({lineage_of[g] for o in oofs.values() for g in o["genome_id"]}, key=str)
    index = {lin: i for i, lin in enumerate(lineages)}
    prepared = {}
    for arm, o in oofs.items():
        prepared[arm] = [(d["y"].to_numpy(), d["p"].to_numpy(),
                          np.array([index[lineage_of[g]] for g in d["genome_id"]]))
                         for _, d in o.groupby("repeat", sort=True)]
    rng = np.random.default_rng(seed)
    draws = {arm: np.empty(n_boot) for arm in oofs}
    for b in range(n_boot):
        counts = np.bincount(rng.integers(0, len(lineages), len(lineages)),
                             minlength=len(lineages))
        for arm, reps in prepared.items():
            draws[arm][b] = np.nanmean([_auc(y, p, counts[c]) for y, p, c in reps])
    lo, hi = (1 - level) / 2 * 100, (1 + level) / 2 * 100

    def interval(x):
        x = x[~np.isnan(x)]
        return {"low": float(np.percentile(x, lo)), "high": float(np.percentile(x, hi)),
                "n_valid": int(x.size)}

    out = {arm: interval(v) for arm, v in draws.items()}
    if diff is not None:
        out[f"{diff[0]}_minus_{diff[1]}"] = interval(draws[diff[0]] - draws[diff[1]])
    out.update({"n_boot": n_boot, "seed": seed, "level": level, "n_lineages": len(lineages)})
    return out
