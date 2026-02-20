import argparse
import pandas as pd
from utils import read_yaml

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--q1q4", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    pairs = cfg["paper6"]["competing_pairs"]

    q = pd.read_csv(args.q1q4)

    top10 = q.sort_values("delta_Q4_minus_Q1", ascending=False).head(10).copy()
    top10["rank_delta"] = range(1, top10.shape[0] + 1)

    want = [f"{p['ligand']}→{p['receptor']}" for p in pairs]
    dom = q[q["pair"].isin(want)].copy()
    dom = dom.sort_values("delta_Q4_minus_Q1", ascending=False)
    dom["rank_in_all_by_delta"] = dom["delta_Q4_minus_Q1"].rank(ascending=False, method="min").astype(int)

    out = pd.concat([
        pd.DataFrame([{"section": "TOP10_Q4_MINUS_Q1"}]),
        top10,
        pd.DataFrame([{"section": "COMPETING_PAIRS"}]),
        dom
    ], ignore_index=True)

    out.to_csv(args.out, index=False)
    print(f"[wrote] {args.out}")

if __name__ == "__main__":
    main()
