import os
import sys
import yaml
import numpy as np
import pandas as pd

def read_yaml(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def die(msg: str, code: int = 2) -> None:
    print(f"[error] {msg}", file=sys.stderr)
    raise SystemExit(code)

def get_expr_vec(ad, gene: str, use_raw: bool = True, layer: str | None = None) -> np.ndarray:
    """Return dense 1D expression vector for gene in AnnData.
    Prefers .raw if use_raw and available, else uses specified layer, else .X.
    """
    import scipy.sparse as sp

    if use_raw and ad.raw is not None and gene in ad.raw.var_names:
        X = ad.raw[:, gene].X
    else:
        if gene not in ad.var_names:
            die(f"Gene {gene} not found in ad.var_names (and not in ad.raw if enabled).")
        if layer is not None:
            if layer not in ad.layers:
                die(f"Requested layer '{layer}' not in ad.layers.")
            X = ad.layers[layer][:, ad.var_names.get_loc(gene)]
        else:
            X = ad[:, gene].X

    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X).reshape(-1)
    return X

def spearman_perm(x: np.ndarray, y: np.ndarray, n_perm: int = 1000, seed: int = 0) -> tuple[float, float]:
    """Spearman correlation with permutation p-value (two-sided)."""
    from scipy.stats import spearmanr
    rng = np.random.default_rng(seed)

    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if x.size < 3:
        return np.nan, np.nan

    r_obs, _ = spearmanr(x, y)
    if np.isnan(r_obs):
        return np.nan, np.nan

    r_perm = []
    y0 = y.copy()
    for _ in range(n_perm):
        rng.shuffle(y0)
        r, _ = spearmanr(x, y0)
        if not np.isnan(r):
            r_perm.append(r)
    r_perm = np.asarray(r_perm)
    if r_perm.size == 0:
        return float(r_obs), np.nan

    p = (np.sum(np.abs(r_perm) >= np.abs(r_obs)) + 1) / (r_perm.size + 1)
    return float(r_obs), float(p)

def qcut_quartiles(s: pd.Series) -> pd.Series:
    """Robust quartiles (labels 1..4). Falls back to rank-based if ties cause issues."""
    try:
        return pd.qcut(s, 4, labels=[1, 2, 3, 4])
    except Exception:
        r = s.rank(method="average")
        return pd.qcut(r, 4, labels=[1, 2, 3, 4])
