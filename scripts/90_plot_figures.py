import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utils import ensure_dir, read_yaml

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--all", required=True)
    ap.add_argument("--q1q4", required=True)
    ap.add_argument("--wilc", required=True)
    ap.add_argument("--corr", required=True)
    ap.add_argument("--spatial", required=True)
    ap.add_argument("--dom", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    _ = read_yaml(args.config)
    ensure_dir(args.outdir)

    q1q4 = pd.read_csv(args.q1q4)
    dom = pd.read_csv(args.dom)
    corr = pd.read_csv(args.corr)
    spatial = pd.read_csv(args.spatial)
    wilc = pd.read_csv(args.wilc)

    top = q1q4.sort_values("delta_Q4_minus_Q1", ascending=False).head(15)
    plt.figure()
    plt.barh(top["pair"][::-1], top["delta_Q4_minus_Q1"][::-1])
    plt.xlabel("Delta strength (Q4 - Q1)")
    plt.title("Activated vs resting B: enriched B→T interactions")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "bar_top15_q4_minus_q1.png"), dpi=200)
    plt.close()

    plt.figure()
    plt.axis("off")
    txt = (
        f"{corr['pair'].iloc[0]}\n"
        f"n={int(corr['n_samples'].iloc[0])}\n"
        f"Spearman r Th2={corr['spearman_r_Th2'].iloc[0]:.3f}, p={corr['perm_p_Th2'].iloc[0]:.3g}\n"
        f"Spearman r Th17={corr['spearman_r_Th17'].iloc[0]:.3f}, p={corr['perm_p_Th17'].iloc[0]:.3g}"
    )
    plt.text(0.02, 0.98, txt, va="top")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "text_corr_summary.png"), dpi=200)
    plt.close()

    comp = dom[dom.get("section").isna() & dom["pair"].notna()].copy()
    if not comp.empty and "delta_Q4_minus_Q1" in comp.columns:
        plt.figure()
        plt.barh(comp["pair"][::-1], comp["delta_Q4_minus_Q1"][::-1])
        plt.xlabel("Delta strength (Q4 - Q1)")
        plt.title("Competing suppressive interactions (activated B→T)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "bar_competing_pairs.png"), dpi=200)
        plt.close()

    plt.figure()
    plt.axis("off")
    s = spatial.iloc[0].to_dict()
    txt2 = (
        f"Immune-enriched spots: {int(s['n_spots_immune'])}\n"
        f"Spearman(LGALS9,HAVCR2) r={s['spearman_r_LGALS9_vs_HAVCR2']:.3f}, p={s['perm_p_LGALS9_vs_HAVCR2']:.3g}\n"
        f"LGALS9-high spots: {int(s['n_LGALS9_high_spots'])}\n"
        f"HAVCR2-high spots: {int(s['n_HAVCR2_high_spots'])}\n"
        f"Mean dist LGALS9hi→nearest HAVCR2hi: {s['mean_dist_LGALS9hi_to_nearest_HAVCR2hi']:.2f}\n"
        f"Median dist LGALS9hi→nearest HAVCR2hi: {s['median_dist_LGALS9hi_to_nearest_HAVCR2hi']:.2f}"
    )
    plt.text(0.02, 0.98, txt2, va="top")
    plt.tight_layout()
    plt.savefig(os.path.join(args.outdir, "text_spatial_summary.png"), dpi=200)
    plt.close()

    w_summary = wilc[wilc["row_type"] == "summary"].iloc[0]
    w_pairs = wilc[wilc["row_type"] == "per_sample"].copy()

    if not w_pairs.empty:
        plt.figure()
        for _, row in w_pairs.iterrows():
            plt.plot([0, 1], [row["strength_Q1"], row["strength_Q4"]], marker="o")
        plt.xticks([0, 1], ["Q1 (resting)", "Q4 (activated)"])
        plt.ylabel("Interaction strength (−aggregate_rank)")
        plt.title(
            f"{w_summary['pair']} paired Wilcoxon (Q4>Q1): "
            f"p={w_summary['wilcoxon_p_greater_Q4']:.3g}, n={int(w_summary['n_samples_paired'])}"
        )
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "lgals9_havcr2_q4_vs_q1_wilcoxon.jpg"), dpi=300)
        plt.close()

    print(f"[ok] figures written to {args.outdir}")

if __name__ == "__main__":
    main()
