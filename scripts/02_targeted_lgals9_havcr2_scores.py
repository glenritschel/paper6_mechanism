#!/usr/bin/env python
import argparse
import numpy as np
import pandas as pd
import scanpy as sc

from utils import read_yaml, die, get_expr_vec


def q1_q4_positions_by_rank(vals: np.ndarray, positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Given values (same length as positions), select bottom/top 25% by rank.
    Robust to ties. Returns (q1_positions, q4_positions) as integer indices into AnnData.
    """
    v = pd.to_numeric(pd.Series(vals), errors="coerce").to_numpy()
    ok = np.isfinite(v)
    v = v[ok]
    pos = positions[ok]

    n = v.shape[0]
    if n < 8:
        return (np.array([], dtype=int), np.array([], dtype=int))

    # rank (average ties), then pick bottom/top k
    r = pd.Series(v).rank(method="average").to_numpy()
    k = max(1, int(np.floor(0.25 * n)))

    q1 = pos[np.argsort(r)[:k]]
    q4 = pos[np.argsort(r)[-k:]]
    return (q1.astype(int), q4.astype(int))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]
    tcfg = c.get("targeted", {})

    ad = sc.read_h5ad(c["adata_path"])

    ct = c["cell_type_key"]
    sk = c["sample_key"]
    b_label = c["b_cell_label"]
    t_label = c["t_cell_label"]
    b_act = c["b_activation_key"]

    ligand = tcfg.get("ligand_gene", "LGALS9")
    receptor = tcfg.get("receptor_gene", "HAVCR2")

    min_bq = int(tcfg.get("min_b_cells_per_quartile", 5))
    min_t = int(tcfg.get("min_t_cells", 10))

    if ligand not in ad.var_names:
        die(f"Targeted ligand gene not in var_names: {ligand}")
    if receptor not in ad.var_names:
        die(f"Targeted receptor gene not in var_names: {receptor}")
    if b_act not in ad.obs.columns:
        die(f"Missing B activation obs key: {b_act}")

    # Expression vectors across ALL cells (dense 1D)
    lg_vec = get_expr_vec(ad, ligand, use_raw=bool(c.get("use_raw", False)), layer=c.get("layer"))
    rc_vec = get_expr_vec(ad, receptor, use_raw=bool(c.get("use_raw", False)), layer=c.get("layer"))

    out_rows = []

    # group -> integer positions
    for sid, idx in ad.obs.groupby(sk, sort=False).indices.items():
        idx = np.asarray(idx, dtype=int)

        ct_s = ad.obs.iloc[idx][ct].astype(str).to_numpy()
        bmask = (ct_s == b_label)
        tmask = (ct_s == t_label)

        t_pos = idx[tmask]
        if t_pos.size < min_t:
            continue

        b_pos = idx[bmask]
        if b_pos.size < (min_bq * 2):
            continue

        b_act_vals = pd.to_numeric(ad.obs.iloc[b_pos][b_act], errors="coerce").to_numpy()
        bq1_pos, bq4_pos = q1_q4_positions_by_rank(b_act_vals, b_pos)

        if bq1_pos.size < min_bq or bq4_pos.size < min_bq:
            continue

        lg_q1 = float(np.mean(lg_vec[bq1_pos]))
        lg_q4 = float(np.mean(lg_vec[bq4_pos]))
        rc_t = float(np.mean(rc_vec[t_pos]))

        out_rows.append({
            "sample_id": sid,
            "n_B_Q1": int(bq1_pos.size),
            "n_B_Q4": int(bq4_pos.size),
            "n_T": int(t_pos.size),
            "mean_LGALS9_B_Q1": lg_q1,
            "mean_LGALS9_B_Q4": lg_q4,
            "mean_HAVCR2_T": rc_t,
            "Q1_strength": lg_q1 * rc_t,
            "Q4_strength": lg_q4 * rc_t,
        })

    df = pd.DataFrame(out_rows)
    if df.empty:
        die("No samples passed targeted scoring filters. Lower min thresholds or confirm labels/keys.")
    df = df.sort_values("sample_id")
    df.to_csv(args.out, index=False)
    print(f"[wrote] {args.out} rows={df.shape[0]}")


if __name__ == "__main__":
    main()

