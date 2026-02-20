import argparse
import pandas as pd
import scanpy as sc
from utils import read_yaml, spearman_perm, die

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--liana", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nperm", type=int, default=1000)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]
    pair = cfg["paper6"]["pair_name"]

    ad = sc.read_h5ad(c["adata_path"])
    sample_key = c["sample_key"]
    th2_key = c["th2_key"]
    th17_key = c["th17_key"]
    ct = c["cell_type_key"]
    t_label = c["t_cell_label"]

    tmask = (ad.obs[ct].astype(str) == t_label)
    prog = ad.obs.loc[tmask, [sample_key, th2_key, th17_key]].copy()
    prog[th2_key] = pd.to_numeric(prog[th2_key], errors="coerce")
    prog[th17_key] = pd.to_numeric(prog[th17_key], errors="coerce")
    prog_agg = prog.groupby(sample_key).median(numeric_only=True).reset_index()

    lr = pd.read_csv(args.liana)
    if "pair" not in lr.columns:
        die("LIANA table is missing 'pair' column. Re-run Step 1.")

    lower = {x.lower(): x for x in lr.columns}
    sample_col = None
    for cand in ["sample", "sample_id", "patient", "donor"]:
        if cand in lower:
            sample_col = lower[cand]
            break

    if "aggregate_rank" not in lower:
        die("Expected 'aggregate_rank' in LIANA output for Paper 4 consistency.")
    ar_col = lower["aggregate_rank"]

    if sample_col is None:
        die("LIANA output does not have a sample column; cannot do per-sample correlation. "
            "Ensure Step 1 outputs per-sample LIANA rows.")

    sub = lr[lr["pair"] == pair].copy()
    if sub.empty:
        die(f"Pair {pair} not found in LIANA table. Check naming (complex vs gene).")

    sub["interaction_strength"] = -pd.to_numeric(sub[ar_col], errors="coerce")
    sub = sub[[sample_col, "interaction_strength"]].rename(columns={sample_col: sample_key})

    df = prog_agg.merge(sub, on=sample_key, how="inner")
    if df.shape[0] < 5:
        die(f"Too few samples with both program scores and interaction strength: {df.shape[0]}")

    r_th2, p_th2 = spearman_perm(df["interaction_strength"].to_numpy(), df[th2_key].to_numpy(),
                                n_perm=args.nperm, seed=0)
    r_th17, p_th17 = spearman_perm(df["interaction_strength"].to_numpy(), df[th17_key].to_numpy(),
                                  n_perm=args.nperm, seed=1)

    out = pd.DataFrame([{
        "pair": pair,
        "n_samples": int(df.shape[0]),
        "spearman_r_Th2": r_th2,
        "perm_p_Th2": p_th2,
        "spearman_r_Th17": r_th17,
        "perm_p_Th17": p_th17
    }])
    out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

if __name__ == "__main__":
    main()
