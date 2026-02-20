import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from utils import read_yaml, get_expr_vec, spearman_perm, die

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nperm", type=int, default=1000)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse249279"]
    ligand = cfg["paper6"]["sender_gene"]
    receptor = cfg["paper6"]["receiver_gene"]

    ad = sc.read_h5ad(c["adata_path"])
    use_raw = bool(c.get("use_raw", True))
    layer = c.get("layer", None)

    f = c["immune_filter"]
    if f["mode"] == "quantile":
        k = f["obs_key"]

        if k not in ad.obs.columns:
            # allow a computed proxy based on existing score columns
            if k == "immune_proxy":
                needed = ["T_score", "B_score", "Plasma_score"]
                missing = [x for x in needed if x not in ad.obs.columns]
                if missing:
                    die(f"Cannot compute immune_proxy; missing obs columns: {missing}")
                ad.obs["immune_proxy"] = (
                    ad.obs["T_score"].astype(float)
                    + ad.obs["B_score"].astype(float)
                    + ad.obs["Plasma_score"].astype(float)
                )
            else:
                die(f"Spatial immune score obs key missing: {k}")

        q = float(f["quantile"])
        vals = pd.to_numeric(ad.obs[k], errors="coerce").to_numpy()
        thr = np.nanquantile(vals, q)
        keep = vals >= thr
    elif f["mode"] == "obs_bool":
        k = f["obs_bool_key"]
        if k not in ad.obs.columns:
            die(f"Spatial immune bool obs key missing: {k}")
        keep = ad.obs[k].astype(bool).to_numpy()
    else:
        die(f"Unknown immune_filter.mode: {f['mode']}")

    ad2 = ad[keep].copy()
    if ad2.n_obs < 10:
        die(f"Too few immune-enriched spots after filter: {ad2.n_obs}")

    x = get_expr_vec(ad2, ligand, use_raw=use_raw, layer=layer)
    y = get_expr_vec(ad2, receptor, use_raw=use_raw, layer=layer)

    r, p = spearman_perm(x, y, n_perm=args.nperm, seed=0)

    coords_key = c.get("coords_key", "spatial")
    if coords_key not in ad2.obsm:
        die(f"Missing spatial coordinates in ad.obsm['{coords_key}']")
    coords = np.asarray(ad2.obsm[coords_key])

    qexpr = float(c.get("top_quantile_expr", 0.90))
    x_thr = np.nanquantile(x, qexpr)
    y_thr = np.nanquantile(y, qexpr)

    x_hi = coords[x >= x_thr]
    y_hi = coords[y >= y_thr]

    def nearest_dist(A, B):
        if A.shape[0] == 0 or B.shape[0] == 0:
            return np.array([])
        d = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(axis=2))
        return d.min(axis=1)

    dists = nearest_dist(x_hi, y_hi)

    out = pd.DataFrame([{
        "n_spots_immune": int(ad2.n_obs),
        "spearman_r_LGALS9_vs_HAVCR2": r,
        "perm_p_LGALS9_vs_HAVCR2": p,
        "q_expr": qexpr,
        "n_LGALS9_high_spots": int(x_hi.shape[0]),
        "n_HAVCR2_high_spots": int(y_hi.shape[0]),
        "mean_dist_LGALS9hi_to_nearest_HAVCR2hi": float(np.nanmean(dists)) if dists.size else np.nan,
        "median_dist_LGALS9hi_to_nearest_HAVCR2hi": float(np.nanmedian(dists)) if dists.size else np.nan
    }])
    out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

if __name__ == "__main__":
    main()
