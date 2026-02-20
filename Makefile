SHELL := /bin/bash
PY := python
CFG := configs/paper6.yaml

REPORTS := reports
FIGS := reports/figures

.PHONY: all check clean

all: \
  $(REPORTS)/liana_bt_interactions_all.csv \
  $(REPORTS)/liana_q1_vs_q4_comparison.csv \
  $(REPORTS)/lgals9_havcr2_q4_vs_q1_wilcoxon.csv \
  $(REPORTS)/lgals9_havcr2_vs_Th2_Th17_correlation.csv \
  $(REPORTS)/spatial_lgals9_havcr2_colocalization.csv \
  $(REPORTS)/top10_competing_interactions.csv \
  $(FIGS)/DONE

check:
	$(PY) scripts/00_check_inputs.py --config $(CFG)

$(REPORTS):
	mkdir -p $(REPORTS)

$(FIGS):
	mkdir -p $(FIGS)

$(REPORTS)/liana_bt_interactions_all.csv: scripts/01_liana_bt_all.py scripts/utils.py $(CFG) | $(REPORTS)
	$(PY) $< --config $(CFG) --out $@

$(REPORTS)/liana_q1_vs_q4_comparison.csv: scripts/02_liana_q1_q4.py scripts/utils.py $(CFG) | $(REPORTS)
	$(PY) $< --config $(CFG) --out $@

$(REPORTS)/lgals9_havcr2_q4_vs_q1_wilcoxon.csv: scripts/02b_wilcoxon_q1_q4_pair.py scripts/utils.py $(CFG) | $(REPORTS)
	$(PY) $< --config $(CFG) --out $@

$(REPORTS)/lgals9_havcr2_vs_Th2_Th17_correlation.csv: scripts/03_corr_lgals9_tim3_vs_programs.py scripts/utils.py $(CFG) \
  $(REPORTS)/liana_bt_interactions_all.csv | $(REPORTS)
	$(PY) $< --config $(CFG) --liana $(REPORTS)/liana_bt_interactions_all.csv --out $@

$(REPORTS)/spatial_lgals9_havcr2_colocalization.csv: scripts/04_spatial_lgals9_havcr2.py scripts/utils.py $(CFG) | $(REPORTS)
	$(PY) $< --config $(CFG) --out $@

$(REPORTS)/top10_competing_interactions.csv: scripts/05_competing_interactions.py scripts/utils.py $(CFG) \
  $(REPORTS)/liana_q1_vs_q4_comparison.csv | $(REPORTS)
	$(PY) $< --config $(CFG) --q1q4 $(REPORTS)/liana_q1_vs_q4_comparison.csv --out $@

$(FIGS)/DONE: scripts/90_plot_figures.py scripts/utils.py $(CFG) \
  $(REPORTS)/liana_bt_interactions_all.csv \
  $(REPORTS)/liana_q1_vs_q4_comparison.csv \
  $(REPORTS)/lgals9_havcr2_q4_vs_q1_wilcoxon.csv \
  $(REPORTS)/lgals9_havcr2_vs_Th2_Th17_correlation.csv \
  $(REPORTS)/spatial_lgals9_havcr2_colocalization.csv \
  $(REPORTS)/top10_competing_interactions.csv | $(FIGS)
	$(PY) $< --config $(CFG) \
	  --all $(REPORTS)/liana_bt_interactions_all.csv \
	  --q1q4 $(REPORTS)/liana_q1_vs_q4_comparison.csv \
	  --wilc $(REPORTS)/lgals9_havcr2_q4_vs_q1_wilcoxon.csv \
	  --corr $(REPORTS)/lgals9_havcr2_vs_Th2_Th17_correlation.csv \
	  --spatial $(REPORTS)/spatial_lgals9_havcr2_colocalization.csv \
	  --dom $(REPORTS)/top10_competing_interactions.csv \
	  --outdir $(FIGS)
	touch $@

clean:
	rm -rf reports/*
