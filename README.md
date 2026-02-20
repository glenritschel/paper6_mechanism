# Paper 6 mechanism pipeline (SSc)

Makefile-driven, idempotent pipeline for:

1) LIANA B→T interactions (aggregate_rank)
2) Q1 vs Q4 B_activation stratification
3) Paired Wilcoxon test for LGALS9→HAVCR2 (Q4 > Q1)
4) Sample-level correlation (LGALS9→HAVCR2 strength vs Th2/Th17)
5) Spatial cross-reference in GSE249279 (immune-enriched spots)
6) Figures

## Quick start

1. Edit `configs/paper6.yaml` to point to your processed AnnData paths and correct obs keys.
2. Run:

```bash
make check
make -j 4
```

Outputs appear in `reports/` and `reports/figures/`.
