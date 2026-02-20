import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from utils import read_yaml, qcut_quartiles, die

def run_liana(ad, groupby: str, use_raw: bool):
    import liana as li
    li.mt.rank_aggregate(ad, groupby=groupby, use_raw=use_raw, verbose=False)
    for k in ["liana_res", "liana_results", "rank_aggregate"]:
        if k in ad.uns:
            return pd.DataFrame(ad.uns[k])
    die("Could not find LIANA results in ad.uns after rank_aggregate.")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]
    pair = cfg["paper6"]["pair_name"]

    ad = sc.read_h5ad(c["adata_path"])

    ct = c["cell_type_key"]
    sample_key = c["sample_key"]
    b_act = c["b_activation_key"]
    b_label = c["b_cell_label"]
    t_label = c["t_cell_label"]
    use_raw = bool(c.get("use_raw", True))

    bmask = (ad.obs[ct].astype(str) == b_label)
    if bmask.sum() == 0:
        die("No B cells found.")

    q = qcut_quartiles(pd.to_numeric(ad.obs.loc[bmask, b_act], errors="coerce"))
    ad.obs.loc[bmask, "B_quartile"] = q.astype("Int64")

    keep_bt = ad.obs[ct].astype(str).isin([b_label, t_label])
    ad_bt = ad[keep_bt].copy()

    min_s = int(c.get("min_cells_per_sample_sender", 3))
    min_r = int(c.get("min_cells_per_sample_receiver", 5))

    bq = ad_bt.obs.loc[ad_bt.obs[ct].astype(str) == b_label, [sample_key, "B_quartile"]].copy()
    tcounts = ad_bt.obs.loc[ad_bt.obs[ct].astype(str) == t_label].groupby(sample_key).size()

    bq_counts = bq.groupby([sample_key, "B_quartile"]).size().unstack(fill_value=0)
    eligible = []
    for sid, row in bq_counts.iterrows():
        if row.get(1, 0) >= min_s and row.get(4, 0) >= min_s and tcounts.get(sid, 0) >= min_r:
            eligible.append(sid)

    if len(eligible) < 5:
        die(f"Too few eligible samples for paired test: {len(eligible)}")

    rows = []
    for sid in eligible:
        ad_s = ad_bt[ad_bt.obs[sample_key] == sid].copy()

        ad_s.obs["paper6_group"] = ad_s.obs[ct].astype(str)
        bmask_s = (ad_s.obs[ct].astype(str) == b_label)

        q_s = ad_s.obs["B_quartile"]
        ad_s.obs.loc[bmask_s & (q_s == 1), "paper6_group"] = "B_Q1"
        ad_s.obs.loc[bmask_s & (q_s == 4), "paper6_group"] = "B_Q4"

        keep = ad_s.obs["paper6_group"].isin(["B_Q1", "B_Q4", t_label])
        ad_s = ad_s[keep].copy()

        # Guard: ensure each required group has >= min cells AFTER we subset to B_Q1/B_Q4/T
        counts = ad_s.obs["paper6_group"].value_counts()
        n_q1 = int(counts.get("B_Q1", 0))
        n_q4 = int(counts.get("B_Q4", 0))
        n_t  = int(counts.get(t_label, 0))

        if n_q1 < min_s or n_q4 < min_s or n_t < min_r:
            continue

        # Additional guard: if any group is completely absent, LIANA logFC can blow up
        if n_q1 == 0 or n_q4 == 0 or n_t == 0:
            continue

        try:
            res = run_liana(ad_s, groupby="paper6_group", use_raw=use_raw)
        except ZeroDivisionError:
            # LIANA can hit zero-denom in logFC when "rest" group is empty
            continue
        except Exception:
            # If any other per-sample LIANA failure occurs, skip this sample
            continue

        if "aggregate_rank" not in res.columns:
            if ("magnitude_rank" in res.columns) and ("specificity_rank" in res.columns):
                res["aggregate_rank"] = (pd.to_numeric(res["magnitude_rank"], errors="coerce")
                                       + pd.to_numeric(res["specificity_rank"], errors="coerce")) / 2.0
            elif "lrscore" in res.columns:
                res["aggregate_rank"] = pd.to_numeric(res["lrscore"], errors="coerce").rank(ascending=False, method="average")
            else:
                die(f"Cannot construct aggregate_rank in per-sample run; columns: {list(res.columns)}")

        cols = {x.lower(): x for x in res.columns}
        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        ligand_c = col("ligand", "ligand_complex")
        receptor_c = col("receptor", "receptor_complex")
        source_c  = col("source", "source_cell_type", "source_group")
        target_c  = col("target", "target_cell_type", "target_group")
        rank_c    = col("aggregate_rank")

        if not all([ligand_c, receptor_c, source_c, target_c, rank_c]):
            die(f"Unexpected LIANA columns in per-sample run. Columns: {list(res.columns)}")

        res["pair"] = res[ligand_c].astype(str) + "→" + res[receptor_c].astype(str)

        def strength(src):
            hit = res[(res[source_c].astype(str) == src) &
                      (res[target_c].astype(str) == t_label) &
                      (res["pair"] == pair)]
            if hit.shape[0] == 0:
                return np.nan
            ar = pd.to_numeric(hit[rank_c].iloc[0], errors="coerce")
            return -float(ar) if pd.notna(ar) else np.nan

        s_q1 = strength("B_Q1")
        s_q4 = strength("B_Q4")
        if np.isnan(s_q1) or np.isnan(s_q4):
            continue

        rows.append({
            sample_key: sid,
            "pair": pair,
            "strength_Q1": s_q1,
            "strength_Q4": s_q4,
            "delta_Q4_minus_Q1": s_q4 - s_q1
        })

    df = pd.DataFrame(rows)
    if df.shape[0] < 5:
        die(f"Too few paired samples after extracting pair scores: {df.shape[0]}")

    from scipy.stats import wilcoxon
    stat, p = wilcoxon(df["strength_Q4"], df["strength_Q1"], alternative="greater", zero_method="wilcox")

    summary = pd.DataFrame([{
        "pair": pair,
        "n_samples_paired": int(df.shape[0]),
        "wilcoxon_stat": float(stat),
        "wilcoxon_p_greater_Q4": float(p),
        "median_strength_Q1": float(np.median(df["strength_Q1"])),
        "median_strength_Q4": float(np.median(df["strength_Q4"])),
        "median_delta_Q4_minus_Q1": float(np.median(df["delta_Q4_minus_Q1"]))
    }])

    df_out = pd.concat(
        [summary.assign(row_type="summary"),
         df.assign(row_type="per_sample")],
        ignore_index=True
    )

    df_out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out} rows={df_out.shape[0]}")

if __name__ == "__main__":
    main()
