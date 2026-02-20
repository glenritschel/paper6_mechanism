import os
import scanpy as sc
import anndata as ad

# INPUTS (adjust if needed)
base = "../paper_05_bt_cells_spatial/data/processed/spatial_scored"
sections = ["A", "B", "C", "D"]

# OUTPUT (NO tilde; expand explicitly)
outdir = os.path.expanduser("~/ssc_data/processed")
os.makedirs(outdir, exist_ok=True)
outpath = os.path.join(outdir, "GSE249279_spatial_merged.h5ad")

ads = []
for s in sections:
    p = os.path.join(base, f"{s}.h5ad")
    a = sc.read_h5ad(p)

    # track section
    a.obs["section"] = s

    # make obs_names unique across sections (barcode collisions are normal)
    a.obs_names = [f"{s}_{idx}" for idx in a.obs_names]

    ads.append(a)

merged = ad.concat(ads, join="outer", label="section", fill_value=0, index_unique=None)

# one more safety pass
merged.obs_names_make_unique()

merged.write_h5ad(outpath)
print(f"[wrote] {outpath}  n_obs={merged.n_obs}  n_vars={merged.n_vars}")
