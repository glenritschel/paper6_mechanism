import argparse
import pandas as pd
import scanpy as sc
from utils import read_yaml, die


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
    groupby = c["cell_type_key"]
    sample_key = c["sample_key"]
    use_raw = bool(c.get("use_raw", True))

    res = run_liana(ad, groupby=groupby, use_raw=use_raw)

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
    score_c   = col("aggregate_rank", "lrscore", "magnitude_rank", "score")

    if not all([ligand_c, receptor_c, source_c, target_c, score_c]):
        die(f"Unexpected LIANA columns. Found: {list(res.columns)}")

    b = c["b_cell_label"]
    t = c["t_cell_label"]
    bt = res[(res[source_c].astype(str) == b) & (res[target_c].astype(str) == t)].copy()

    pcol = col("pvalue", "p_val", "p_val_adj", "pval", "p_value", "cellphone_pvals")
    if pcol is not None:
        bt = bt[pd.to_numeric(bt[pcol], errors="coerce") < float(c["liana"]["pval_cutoff"])]

    bt["pair"] = bt[ligand_c].astype(str) + "→" + bt[receptor_c].astype(str)

    if score_c.lower() == "aggregate_rank":
        bt = bt.sort_values(by=score_c, ascending=True)
    else:
        bt = bt.sort_values(by=score_c, ascending=False)

    bt.to_csv(args.out, index=False)
    print(f"[wrote] {args.out} rows={bt.shape[0]}")

if __name__ == "__main__":
    main()
