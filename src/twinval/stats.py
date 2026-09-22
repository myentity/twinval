"""Statistical building blocks shared by all analysis levels."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps
from statsmodels.stats.multitest import multipletests


# ------------------------------------------------------------------ agreement / calibration
def ccc_components(H: np.ndarray, X: np.ndarray) -> dict[str, np.ndarray]:
    """Lin's CCC and its decomposition CCC = r * C_b, vectorised over columns.

    C_b is computed directly from the location shift ``u`` and scale ratio ``v``
    (not as CCC / r), which stays stable when r is close to zero::

        v = SD_X / SD_H,  u = (mean_X - mean_H) / sqrt(SD_H * SD_X),  C_b = 2 / (v + 1/v + u^2)
    """
    H = np.atleast_2d(H.T).T if H.ndim == 1 else H
    X = np.atleast_2d(X.T).T if X.ndim == 1 else X
    mh, mx = H.mean(0), X.mean(0)
    sh, sx = H.std(0, ddof=1), X.std(0, ddof=1)
    cov = ((H - mh) * (X - mx)).sum(0) / (H.shape[0] - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = cov / (sh * sx)
        v = sx / sh
        u = (mx - mh) / np.sqrt(sh * sx)
        cb = 2 / (v + 1 / v + u ** 2)
    return {"r": r, "ccc": r * cb, "cb": cb, "u": u, "v": v, "mean_h": mh, "mean_x": mx, "sd_h": sh, "sd_x": sx}


def ccc_1d(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's CCC of two vectors (ddof = 1)."""
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(ddof=1), y.var(ddof=1)
    cov = ((x - mx) * (y - my)).sum() / (len(x) - 1)
    return float(2 * cov / (vx + vy + (mx - my) ** 2))


def human_z(H: np.ndarray, T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Standardise BOTH sides with the human mean/SD of each scale ((X - mean_H) / SD_H).

    Twin scores are never re-centred on their own mean, so systematic twin bias stays visible.
    """
    mu, sd = H.mean(0), H.std(0, ddof=1)
    return (H - mu) / sd, (T - mu) / sd


# ------------------------------------------------------------------ identification (own vs other twin)
def own_cross_stats(M: np.ndarray, lower_is_better: bool) -> dict:
    """Per-person statistics of an N x N Human_i x Twin_j distance (or similarity) matrix.

    * rank of the own twin (1 = best match); * per-person AUC
    ``P(own better than a random other twin)`` with ties counted 1/2 and the own cell excluded.
    """
    n = M.shape[0]
    own = np.diag(M).copy()
    off = ~np.eye(n, dtype=bool)
    cross_mean = np.array([M[i, off[i]].mean() for i in range(n)])
    order = np.argsort(M if lower_is_better else -M, axis=1)
    rank = np.empty(n, dtype=int)
    for i in range(n):
        rank[i] = np.where(order[i] == i)[0][0] + 1
    better = (own[:, None] < M) if lower_is_better else (own[:, None] > M)
    eq = own[:, None] == M
    np.fill_diagonal(better, False)
    np.fill_diagonal(eq, False)
    auc = (better.sum(1) + 0.5 * eq.sum(1)) / (n - 1)
    return {"own": own, "cross_mean": cross_mean, "rank": rank, "auc": auc,
            "percentile": (n - rank) / (n - 1) * 100, "tie_share": float(eq.sum() / (n * (n - 1)))}


def permutation_auc_p(M: np.ndarray, obs_auc: float, rng: np.random.Generator, n_perm: int,
                      lower_is_better: bool = True) -> tuple[float, float, float]:
    """Label-permutation test of the mean AUC: p = (1 + #{null >= obs}) / (B + 1)."""
    n = M.shape[0]
    null = np.empty(n_perm)
    idx = np.arange(n)
    for b in range(n_perm):
        perm = rng.permutation(n)
        own_p = M[idx, perm]
        better = (own_p[:, None] < M) if lower_is_better else (own_p[:, None] > M)
        eq = own_p[:, None] == M
        mask = np.ones((n, n), dtype=bool)
        mask[idx, perm] = False
        cnt = mask.sum(1)
        null[b] = (((better & mask).sum(1) + 0.5 * (eq & mask).sum(1)) / cnt).mean()
    return float((1 + np.sum(null >= obs_auc)) / (n_perm + 1)), float(null.mean()), float(null.std())


def cohens_dz(diff: np.ndarray) -> float:
    sd = diff.std(ddof=1)
    return float(diff.mean() / sd) if sd > 0 else float("nan")


# ------------------------------------------------------------------ resampling
def boot_indices(rng: np.random.Generator, n_boot: int, n: int) -> np.ndarray:
    return rng.integers(0, n, size=(n_boot, n))


def pct_ci(x: np.ndarray, axis: int = 0, level: float = 95.0) -> tuple[np.ndarray, np.ndarray]:
    a = (100 - level) / 2
    return np.nanpercentile(x, a, axis=axis), np.nanpercentile(x, 100 - a, axis=axis)


def grouped_fold_ids(idx: np.ndarray, n_folds: int, rng: np.random.Generator) -> np.ndarray:
    """Fold assignment for a bootstrap sample: all copies of the same original respondent
    go to the same fold (no train/test leakage of duplicated respondents)."""
    unique_orig, inverse = np.unique(idx, return_inverse=True)
    perm = rng.permutation(len(unique_orig))
    group_fold = np.empty(len(unique_orig), dtype=int)
    group_fold[perm] = np.arange(len(unique_orig)) % n_folds
    return group_fold[inverse]


def bh_adjust(p: np.ndarray | pd.Series) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def bootstrap_p(boots: np.ndarray) -> tuple[float, bool]:
    """Two-sided bootstrap p-value with the (b + 1) / (n + 1) correction.

    ``b`` is the number of replicates on the smaller side of zero. A finite number of replicates
    cannot resolve p = 0, so the smallest reportable value is 2 / (n + 1); the flag says when the
    estimate sits at this resolution floor.
    """
    boots = boots[~np.isnan(boots)]
    n = len(boots)
    if n < 100:
        return float("nan"), False
    b = min(int((boots <= 0).sum()), int((boots >= 0).sum()))
    return float(min(2 * (b + 1) / (n + 1), 1.0)), bool(b == 0)


def cluster_bootstrap_spearman(df: pd.DataFrame, xcol: str, ycol: str, cluster_col: str,
                               n_boot: int, rng: np.random.Generator, min_clusters: int = 4) -> dict | None:
    """Spearman rho with a cluster (e.g. per-instrument) bootstrap CI and bootstrap p-value.

    Scales of one instrument are not independent observations, so whole instruments are resampled.
    With fewer than ``min_clusters`` clusters the bootstrap is unstable; CI and p are then NaN and the
    naive Spearman p-value is kept for reference only.
    """
    sub = df.dropna(subset=[xcol, ycol])
    if len(sub) < 5:
        return None
    rho, p_naive = sps.spearmanr(sub[xcol], sub[ycol])
    clusters = sub[cluster_col].unique()
    out = {"rho": float(rho), "p_naive": float(p_naive), "n": len(sub), "n_clusters": len(clusters),
           "ci_lo": np.nan, "ci_hi": np.nan, "p_boot": np.nan, "p_below_resolution": False, "n_valid_boot": 0}
    if len(clusters) < min_clusters:
        return out
    groups = {c: sub[sub[cluster_col] == c] for c in clusters}
    boots = []
    for _ in range(n_boot):
        draw = rng.choice(clusters, size=len(clusters), replace=True)
        cat = pd.concat([groups[c] for c in draw], ignore_index=True)
        if cat[xcol].nunique() < 2 or cat[ycol].nunique() < 2:
            continue
        boots.append(sps.spearmanr(cat[xcol], cat[ycol])[0])
    boots = np.array(boots)
    out["n_valid_boot"] = len(boots)
    if len(boots) > 100:
        out["ci_lo"], out["ci_hi"] = (float(v) for v in np.percentile(boots, [2.5, 97.5]))
        out["p_boot"], out["p_below_resolution"] = bootstrap_p(boots)
    return out


def fisher_z_mean(r: np.ndarray) -> float:
    """Back-transformed mean of Fisher-z transformed correlations."""
    r = np.asarray(r, float)
    return float(np.tanh(np.arctanh(np.clip(r, -0.999999, 0.999999)).mean()))
