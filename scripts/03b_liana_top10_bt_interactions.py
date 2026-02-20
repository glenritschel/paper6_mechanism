#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--liana_all", required=True, help="CSV from 01_liana_bt_all.py")
    ap.add_argument("--out", required=True, help="Output CSV of topK")
    ap.add_argument("--fig", required=True, help="Output figure path (JPG)")
    ap.add_argument("--topk", type=int, default=10)
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    df = pd.read_csv(args.liana_all)

    # Robust column expectations across LIANA versions
    if "aggregate_rank" not in df.columns:
        raise SystemExit(f"[error] Missing aggregate_rank in LIANA table. Columns: {list(df.columns)}")

    # Pair label
    if "pair" in df.columns:
        df["pair_label"] = df["pair"].astype(str)
    else:
        # LIANA usually uses ligand_complex/receptor_complex
        for c in ["ligand_complex", "receptor_complex"]:
            if c not in df.columns:
                raise SystemExit(f"[error] Missing column in LIANA all table: {c}")
        df["pair_label"] = df["ligand_complex"].astype(str) + "→" + df["receptor_complex"].astype(str)

    # Strength definition (Option A for rank-based)
    # smaller aggregate_rank => stronger; convert to inverted strength for plotting
    df["aggregate_rank"] = pd.to_numeric(df["aggregate_rank"], errors="coerce")
    df = df.dropna(subset=["aggregate_rank"]).copy()
    df["strength"] = 1.0 / (df["aggregate_rank"] + 1e-12)

    top = df.sort_values("aggregate_rank", ascending=True).head(args.topk).copy()
    top = top.sort_values("strength", ascending=True)  # so largest appears at top in barh

    # Write CSV (with explicit strength)
    out_df = top[["pair_label", "aggregate_rank", "strength"]].rename(columns={"pair_label": "pair"})
    out_df.to_csv(args.out, index=False)
    print(f"[wrote] {args.out} rows={out_df.shape[0]}")

    # Plot: horizontal bar plot
    plt.figure(figsize=(8, 4.8))
    plt.barh(out_df["pair"], out_df["strength"])

    # Add numeric labels to bars (optional but helpful)
    for i, v in enumerate(out_df["strength"].values):
        plt.text(v, i, f" {v:.3g}", va="center", fontsize=8)

    plt.xlabel("Interaction strength (1 / aggregate_rank)")
    plt.ylabel("Ligand→Receptor")
    plt.title(f"Top {out_df.shape[0]} B→T interactions (LIANA aggregate_rank)")
    plt.tight_layout()
    plt.savefig(args.fig, dpi=args.dpi)
    plt.close()
    print(f"[wrote] {args.fig}")

if __name__ == "__main__":
    main()

