#!/usr/bin/env python3
import argparse
import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


def get_nested(cfg: dict, paths: list[str], default=None):
    """Return first found nested key among candidate dot-paths."""
    for p in paths:
        cur = cfg
        ok = True
        for k in p.split("."):
            if isinstance(cur, dict) and k in cur:
                cur = cur[k]
            else:
                ok = False
                break
        if ok:
            return cur
    return default


def perm_spearman(x, y, nperm=1000, seed=0):
    rng = np.random.default_rng(seed)
    r_obs, _ = spearmanr(x, y)
    count = 0
    for _ in range(int(nperm)):
        y_perm = rng.permutation(y)
        r_perm, _ = spearmanr(x, y_perm)
        if abs(r_perm) >= abs(r_obs):
            count += 1
    p = (count + 1) / (nperm + 1)
    return float(r_obs), float(p)


def to_dense_1d(X):
    # AnnData slice can be sparse or dense
    if hasattr(X, "toarray"):
        return X.toarray().ravel()
    return np.asarray(X).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fig", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--nperm", type=int, default=1000)
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config))

    # --- Resolve config keys robustly ---
    spatial_path = get_nested(cfg, [
        "spatial.adata_path",
        "gse249279.adata_path",
        "gse249279.merged_adata_path",
    ])
    if not spatial_path:
        raise KeyError("Could not find spatial AnnData path. Expected one of: "
                       "spatial.adata_path, gse249279.adata_path, gse249279.merged_adata_path")

    ligand = get_nested(cfg, [
        "paper6.sender_gene",
        "paper6.ligand",
        "paper6.genes.ligand",
    ], default="LGALS9")

    receptor = get_nested(cfg, [
        "paper6.receiver_gene",
        "paper6.receptor",
        "paper6.genes.receptor",
    ], default="HAVCR2")

    immune_key = get_nested(cfg, [
        "spatial.immune_score_key",
        "gse249279.immune_score_key",
    ], default="immune_score")

    immune_threshold = float(get_nested(cfg, [
        "spatial.immune_threshold",
        "gse249279.immune_threshold",
    ], default=0.0))

    q_expr = float(get_nested(cfg, [
        "spatial.expression_quantile",
        "gse249279.expression_quantile",
        "gse249279.q_expr",
    ], default=0.90))

    # --- Load spatial data ---
    ad = sc.read_h5ad(spatial_path)

    if immune_key not in ad.obs.columns:
        raise KeyError(f"Spatial immune score obs key missing: {immune_key}. "
                       f"Available obs columns: {list(ad.obs.columns)[:50]}")

    if ligand not in ad.var_names or receptor not in ad.var_names:
        raise ValueError(f"Ligand/receptor not in var_names: {ligand}={ligand in ad.var_names}, "
                         f"{receptor}={receptor in ad.var_names}")

    immune_mask = ad.obs[immune_key].astype(float) > immune_threshold
    ad_imm = ad[immune_mask].copy()

    x = to_dense_1d(ad_imm[:, ligand].X)
    y = to_dense_1d(ad_imm[:, receptor].X)

    r, p = perm_spearman(x, y, nperm=args.nperm)

    # quantile thresholds for "high" spots
    x_thr = float(np.quantile(x, q_expr))
    y_thr = float(np.quantile(y, q_expr))

    df = pd.DataFrame({
        "n_spots_immune": [int(ad_imm.n_obs)],
        f"spearman_r_{ligand}_vs_{receptor}": [r],
        f"perm_p_{ligand}_vs_{receptor}": [p],
        "q_expr": [q_expr],
        f"n_{ligand}_high_spots": [int((x >= x_thr).sum())],
        f"n_{receptor}_high_spots": [int((y >= y_thr).sum())],
        f"frac_{ligand}_high": [float((x >= x_thr).mean())],
        f"frac_{receptor}_high": [float((y >= y_thr).mean())],
        "x_quantile_threshold": [x_thr],
        "y_quantile_threshold": [y_thr],
    })

    df.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

    # --- Scatter plot ---
    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, alpha=0.6, s=16)
    plt.axvline(x_thr, linestyle="--")
    plt.axhline(y_thr, linestyle="--")
    plt.xlabel(ligand)
    plt.ylabel(receptor)
    plt.title(f"{ligand} vs {receptor} (immune-enriched spots)\nSpearman r={r:.3f}, perm p={p:.3f}")
    plt.tight_layout()
    plt.savefig(args.fig, dpi=args.dpi)
    plt.close()
    print(f"[wrote] {args.fig}")


if __name__ == "__main__":
    main()
