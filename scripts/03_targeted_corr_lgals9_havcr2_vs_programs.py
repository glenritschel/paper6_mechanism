#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from utils import read_yaml, die, spearman_perm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--nperm", type=int, default=1000)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]

    df = pd.read_csv(args.scores)
    if df.empty:
        die("Scores table is empty.")

    # compute per-sample Th2/Th17 program scores from the PBMC AnnData (T cells only by default)
    ad = sc.read_h5ad(c["adata_path"])
    ct = c["cell_type_key"]
    sk = c["sample_key"]
    t_label = c["t_cell_label"]
    th2_key = c["th2_key"]
    th17_key = c["th17_key"]

    for k in [th2_key, th17_key]:
        if k not in ad.obs.columns:
            die(f"Missing program score obs key: {k}")

    tmask = ad.obs[ct].astype(str) == t_label
    prog = ad.obs.loc[tmask, [sk, th2_key, th17_key]].copy()
    prog[th2_key] = pd.to_numeric(prog[th2_key], errors="coerce")
    prog[th17_key] = pd.to_numeric(prog[th17_key], errors="coerce")

    per_sample = prog.groupby(sk).mean(numeric_only=True).reset_index().rename(columns={
        sk: "sample_id",
        th2_key: "Th2_score",
        th17_key: "Th17_score",
    })

    m = df.merge(per_sample, on="sample_id", how="inner")
    if m.shape[0] < 10:
        die(f"Too few samples after merging scores with programs: {m.shape[0]}")

    x = m["Q4_strength"].astype(float).values
    y_th2 = m["Th2_score"].astype(float).values
    y_th17 = m["Th17_score"].astype(float).values

    r_th2, p_th2 = spearman_perm(x, y_th2, n_perm=args.nperm, seed=0)
    r_th17, p_th17 = spearman_perm(x, y_th17, n_perm=args.nperm, seed=1)

    out = pd.DataFrame([
        {"program": "Th2", "n": int(np.isfinite(x).sum()), "spearman_r": float(r_th2), "perm_p_two_sided": float(p_th2)},
        {"program": "Th17", "n": int(np.isfinite(x).sum()), "spearman_r": float(r_th17), "perm_p_two_sided": float(p_th17)},
    ])
    out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

    # figures
    import os
    import matplotlib.pyplot as plt
    os.makedirs(args.figdir, exist_ok=True)

    def scatter(y, label, r, p, fname):
        plt.figure()
        plt.scatter(x, y, alpha=0.6)
        plt.xlabel("Q4 strength = mean(LGALS9 in B_Q4) × mean(HAVCR2 in T)")
        plt.ylabel(f"{label} program score (mean over T cells)")
        plt.title(f"{label}: Spearman r={r:.3f}, perm p={p:.3g} (n={m.shape[0]})")
        plt.tight_layout()
        plt.savefig(os.path.join(args.figdir, fname), dpi=200)
        plt.close()

    scatter(y_th2, "Th2", r_th2, p_th2, "lgals9_havcr2_vs_th2.jpg")
    scatter(y_th17, "Th17", r_th17, p_th17, "lgals9_havcr2_vs_th17.jpg")
    print(f"[wrote] {args.figdir}/lgals9_havcr2_vs_th2.jpg")
    print(f"[wrote] {args.figdir}/lgals9_havcr2_vs_th17.jpg")


if __name__ == "__main__":
    main()

