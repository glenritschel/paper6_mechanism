import argparse
import pandas as pd
import scanpy as sc
from utils import read_yaml, qcut_quartiles, die


def run_liana(ad, groupby: str, use_raw: bool):
    try:
        import liana as li
    except Exception as e:
        die(f"Failed to import liana. Is LIANA-py installed? Error: {e}")

    # LIANA API varies by version: some accept sample_key, some don't.
    try:
        li.mt.rank_aggregate(
            ad,
            groupby=groupby,
            use_raw=use_raw,
            verbose=True,
        )
    except TypeError:
        # fallback: no sample_key
        li.mt.rank_aggregate(
            ad,
            groupby=groupby,
            use_raw=use_raw,
            verbose=True,
        )

    for k in ["liana_res", "liana_results", "rank_aggregate"]:
        if k in ad.uns:
            return pd.DataFrame(ad.uns[k])
    die("Could not find LIANA results in ad.uns. Check LIANA version / output keys.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]
    ad = sc.read_h5ad(c["adata_path"])

    ct = c["cell_type_key"]
    sample_key = c["sample_key"]
    b_act = c["b_activation_key"]
    b_label = c["b_cell_label"]
    t_label = c["t_cell_label"]
    use_raw = bool(c.get("use_raw", True))

    ad.obs["paper6_group"] = ad.obs[ct].astype(str)

    bmask = (ad.obs[ct].astype(str) == b_label)
    if bmask.sum() == 0:
        die("No B cells found for quartiling.")

    q = qcut_quartiles(pd.to_numeric(ad.obs.loc[bmask, b_act], errors="coerce"))
    ad.obs.loc[bmask, "B_quartile"] = q.astype("Int64")

    ad.obs.loc[bmask & (ad.obs["B_quartile"] == 1), "paper6_group"] = "B_Q1"
    ad.obs.loc[bmask & (ad.obs["B_quartile"] == 4), "paper6_group"] = "B_Q4"

    keep = ad.obs["paper6_group"].isin(["B_Q1", "B_Q4", t_label])
    ad2 = ad[keep].copy()

    min_s = int(c.get("min_cells_per_sample_sender", 20))
    min_r = int(c.get("min_cells_per_sample_receiver", 20))

    vc = ad2.obs.groupby([sample_key, "paper6_group"]).size().unstack(fill_value=0)
    good = []
    for sid, row in vc.iterrows():
        if row.get("B_Q1", 0) >= min_s and row.get("B_Q4", 0) >= min_s and row.get(t_label, 0) >= min_r:
            good.append(sid)
    ad2 = ad2[ad2.obs[sample_key].isin(good)].copy()

    res = run_liana(ad2, groupby="paper6_group", use_raw=use_raw)

    # Ensure aggregate_rank exists for Paper 4 consistency
    if "aggregate_rank" not in res.columns:
        if ("magnitude_rank" in res.columns) and ("specificity_rank" in res.columns):
            res["aggregate_rank"] = (pd.to_numeric(res["magnitude_rank"], errors="coerce")
                                   + pd.to_numeric(res["specificity_rank"], errors="coerce")) / 2.0
        elif "lrscore" in res.columns:
            # fallback: rank lrscore descending (higher lrscore = stronger)
            res["aggregate_rank"] = pd.to_numeric(res["lrscore"], errors="coerce").rank(ascending=False, method="average")
        else:
            die(f"Cannot construct aggregate_rank; columns are: {list(res.columns)}")

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
    score_c   = col("aggregate_rank")
    pcol      = col("pvalue", "p_val", "p_val_adj", "pval", "p_value")

    if not all([ligand_c, receptor_c, source_c, target_c, score_c]):
        die(f"Unexpected LIANA columns. Found: {list(res.columns)}")

    q1 = res[(res[source_c].astype(str) == "B_Q1") & (res[target_c].astype(str) == t_label)].copy()
    q4 = res[(res[source_c].astype(str) == "B_Q4") & (res[target_c].astype(str) == t_label)].copy()

    for df in (q1, q4):
        df["pair"] = df[ligand_c].astype(str) + "→" + df[receptor_c].astype(str)

    m = q4[["pair", score_c] + ([pcol] if pcol else [])].merge(
        q1[["pair", score_c] + ([pcol] if pcol else [])],
        on="pair",
        how="outer",
        suffixes=("_Q4", "_Q1")
    )

    m["strength_Q4"] = -pd.to_numeric(m[f"{score_c}_Q4"], errors="coerce")
    m["strength_Q1"] = -pd.to_numeric(m[f"{score_c}_Q1"], errors="coerce")
    m["delta_Q4_minus_Q1"] = m["strength_Q4"] - m["strength_Q1"]

    if pcol:
        p4 = pd.to_numeric(m.get(f"{pcol}_Q4"), errors="coerce")
        p1 = pd.to_numeric(m.get(f"{pcol}_Q1"), errors="coerce")
        m["sig_either"] = (p4 < float(c["liana"]["pval_cutoff"])) | (p1 < float(c["liana"]["pval_cutoff"]))
    else:
        m["sig_either"] = True

    m = m.sort_values("delta_Q4_minus_Q1", ascending=False)
    m.to_csv(args.out, index=False)
    print(f"[wrote] {args.out} rows={m.shape[0]}")

if __name__ == "__main__":
    main()
