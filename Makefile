SHELL := /bin/bash
PY := python
CFG := configs/paper6.yaml

# ---- CSV outputs ----
LIANA_ALL := reports/liana_bt_interactions_all.csv
LIANA_TOP10 := reports/liana_top10_bt_interactions.csv
TARGETED := reports/targeted_lgals9_havcr2_scores.csv
WILCOXON := reports/lgals9_havcr2_wilcoxon_tim3_detectable.csv
CORR := reports/lgals9_havcr2_vs_Th2_Th17_correlation_tim3_detectable.csv
SPATIAL := reports/spatial_lgals9_havcr2_colocalization.csv

# ---- Figures ----
FIG1 := reports/figures/Figure1_LIANA_Top10_BT_Interactions.jpg
FIG2 := reports/figures/Figure2_LGALS9_HAVCR2_Wilcoxon.jpg
FIG3 := reports/figures/Figure3_Spatial_Colocalization.jpg

.PHONY: all clean dirs figures zip

all: $(LIANA_ALL) $(LIANA_TOP10) $(TARGETED) \
     $(WILCOXON) $(CORR) $(SPATIAL) figures

dirs:
	mkdir -p reports reports/figures

clean:
	rm -rf reports/*

# ---- Step 1: LIANA all ----
$(LIANA_ALL): | dirs
	$(PY) scripts/01_liana_bt_all.py --config $(CFG) --out $@

# ---- Figure 1 + Top10 table ----
$(LIANA_TOP10) $(FIG1): $(LIANA_ALL) | dirs
	$(PY) scripts/03b_liana_top10_bt_interactions.py \
		--liana_all $(LIANA_ALL) \
		--out $(LIANA_TOP10) \
		--fig $(FIG1) \
		--topk 10 --dpi 300

# ---- Targeted scoring ----
$(TARGETED): | dirs
	$(PY) scripts/02_targeted_lgals9_havcr2_scores.py \
		--config $(CFG) --out $@

# ---- Wilcoxon (TIM-3 detectable subset) ----
$(WILCOXON) $(FIG2): $(TARGETED) | dirs
	$(PY) scripts/02b_targeted_lgals9_havcr2_wilcoxon.py \
		--scores $(TARGETED) \
		--out $(WILCOXON) \
		--fig $(FIG2)

# ---- Correlations ----
$(CORR): $(TARGETED) | dirs
	$(PY) scripts/03_targeted_corr_lgals9_havcr2_vs_programs.py \
		--config $(CFG) \
		--scores $(TARGETED) \
		--out $(CORR) \
		--figdir reports/figures \
		--nperm 1000

# ---- Spatial CSV + Figure 3 ----
$(SPATIAL) $(FIG3): | dirs
	$(PY) scripts/04_spatial_lgals9_havcr2.py \
		--config $(CFG) \
		--out $(SPATIAL) \
		--fig $(FIG3) \
		--dpi 300

figures: $(FIG1) $(FIG2) $(FIG3)

zip: figures
	rm -rf paper6_export
	mkdir -p paper6_export
	cp $(FIG1) paper6_export/
	cp $(FIG2) paper6_export/
	cp $(FIG3) paper6_export/
	zip -j paper6_figures.zip paper6_export/*.jpg

