#!/usr/bin/env python
import argparse
import pandas as pd
from scipy.stats import wilcoxon
from utils import die


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True, help="reports/targeted_lgals9_havcr2_scores.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--fig", required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.scores)
    if df.empty:
        die("Scores table is empty.")

    x = df["Q4_strength"].astype(float)
    y = df["Q1_strength"].astype(float)

    keep = x.notna() & y.notna()
    df = df.loc[keep].copy()
    if df.shape[0] < 10:
        die(f"Too few paired samples for Wilcoxon: {df.shape[0]}")

    stat, p = wilcoxon(df["Q4_strength"], df["Q1_strength"], alternative="two-sided", zero_method="wilcox")

    out = pd.DataFrame([{
        "n_pairs": int(df.shape[0]),
        "wilcoxon_stat": float(stat),
        "p_value_two_sided": float(p),
        "median_Q4_strength": float(df["Q4_strength"].median()),
        "median_Q1_strength": float(df["Q1_strength"].median()),
        "median_delta_Q4_minus_Q1": float((df["Q4_strength"] - df["Q1_strength"]).median()),
    }])
    out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

    # figure
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot([0, 1], [df["Q1_strength"].median(), df["Q4_strength"].median()], marker="o")
    # paired points
    for _, r in df.iterrows():
        plt.plot([0, 1], [r["Q1_strength"], r["Q4_strength"]], alpha=0.25)
    plt.xticks([0, 1], ["Q1 (resting)", "Q4 (activated)"])
    plt.ylabel("LGALS9(B) × HAVCR2(T) strength")
    plt.title(f"Paired Wilcoxon p={p:.3g} (n={df.shape[0]})")
    plt.tight_layout()
    plt.savefig(args.fig, dpi=200)
    plt.close()
    print(f"[wrote] {args.fig}")


if __name__ == "__main__":
    main()

