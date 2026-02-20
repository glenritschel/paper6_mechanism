import argparse
import scanpy as sc
from utils import read_yaml, die

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = read_yaml(args.config)
    c = cfg["gse195452"]
    ad = sc.read_h5ad(c["adata_path"])

    required = [c["cell_type_key"], c["sample_key"], c["b_activation_key"], c["th2_key"], c["th17_key"]]
    for k in required:
        if k not in ad.obs.columns:
            die(f"Missing required obs key: {k}")

    ct = c["cell_type_key"]
    labels = set(ad.obs[ct].astype(str))
    if c["b_cell_label"] not in labels:
        die(f"B cell label {c['b_cell_label']} not found in obs[{ct}]")
    if c["t_cell_label"] not in labels:
        die(f"T cell label {c['t_cell_label']} not found in obs[{ct}]")

    print("[ok] GSE195452 AnnData checks passed")

if __name__ == "__main__":
    main()
