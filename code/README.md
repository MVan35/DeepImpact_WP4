# Code for "Neighbouring weed species influence the root microbiota and growth of wheat"

## Getting the data

Download the Zenodo data record (DOI in the article's Data availability statement) and unpack it. Then set, before running:

```
export PAPER_DATA=/path/to/unpacked/data          # folder holding matrice_phase1/, results_v4/, NST/, phylogeny/, traits/
export PAPER_OUT=/path/to/results                 # where the qmd writes (default ../results)
export PAPER_RESULTS=$PAPER_OUT                   # where the notebook reads the qmd outputs
```

Requirements: R ≥ 4.4 with phyloseq, vegan, dplyr, ggplot2, patchwork, car, DESeq2, edgeR, metagenomeSeq, TITAN2, ape, jsonlite;
Python ≥ 3.11 with numpy, pandas, scipy, statsmodels, matplotlib, openpyxl. Runtime about 30 min on a laptop. Versions used for the article: R 4.5 with DESeq2 1.48, edgeR 4.6, metagenomeSeq 1.52. edgeR results change with the edgeR version (4.8.2 gives a larger Fig. 2 call set), so the notebook uses the edgeR exports deposited with the data (`da_exports_paper/`) unless `PAPER_EDGER=regenerated` is set.

Two documents, run in this order from this folder:

1. `Paper_analyses.qmd` (R): the statistical analyses. Writes to `../results/` (stand-alone scripts write to `../results/supporting/`).
2. `Paper_figures.ipynb` (Python, built by `build_notebook.py`): every figure and the Python-side tables. Reads `../results/`,
   writes `figures/`, `tables/`, `differential_abundance/`, and copies each display item to `display_items/` under its
   manuscript name (`Fig_1.png`, `Table_S3.xlsx`, …). Section 3.1.2 calls R for the control-size check.

Both are organised by the Results sections (3.1.1, 3.1.2, 3.2.1, 3.2.2, 3.3); each heading names the display items it produces.

| Display item | Produced in |
|---|---|
| Table 1, Fig. S4; within-soil PERMANOVA (Methods S2) | qmd 3.1.1 |
| Fig. 1, Figs S5, S6, S7 | notebook 3.1.1 (tests: qmd 3.1.1) |
| Fig. 2, Table S3 (with the control-size check) | notebook 3.1.2 (inputs: qmd 3.1.2) |
| Fig. 3, Table S4; Table S5 | notebook 3.2.1; qmd 3.2.1 |
| Fig. 4 (tests: qmd 3.1.1 `revision-analyses`), Table S6 (sensitivity grid), Table S7 (TITAN2) | notebook 3.2.2; qmd 3.2.2 |
| Table 2, Notes S3 | qmd 3.3 (`trait_models_v3.R`) and notebook 3.3 |

Figs S1, S2 and S3 are drawings or images supplied with the manuscript.

## Requirements

* R ≥ 4.2: `phyloseq`, `vegan`, `dplyr`, `ggplot2`, `patchwork`, `car`, `DESeq2`, `edgeR`, `metagenomeSeq`, `TITAN2`, `ape`.
  `V.PhyloMaker2` is optional: the dated tree is supplied (`weeds_wheat_VPM2_WP_S3.tre`) and is rebuilt when the package is installed.
* Python ≥ 3.10: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `openpyxl`.
* Data: `../data`, staged by `stage_submission_data.sh`. The fungal differential-abundance exports and the NST results
  are read from `paper_1/data` until they are staged.

Paths: set `PROJECT` at the top of each document. `PAPER_OUT` (both), `PAPER_RESULTS` and `RSCRIPT` (notebook) redirect
outputs, inputs and the R executable without editing the code.

## Files

| File | Used by | Content |
|---|---|---|
| `validity_checks.R`, `robustness_and_tables.R`, `variance_partition_plot.R` | qmd 3.1.1 | PERMANOVA, assumptions, CLR sensitivity |
| `permanova_order.R` | qmd 3.1.1 | order-independence of Table 1 |
| `revision_analyses.R` | qmd 3.1.1 | within-soil PERMANOVA (Table 1); Fig. 4 db-RDA on species centroids with conditional tests; per-taxon screen; Glomeromycota tests |
| `edger_export.R`, `da_codetection.R` | qmd 3.1.2 | edgeR contrasts against wheat × wheat; export parsing |
| (no module; chunk `table-s5-biomass`) | qmd 3.2.1 | Table S5: log wheat biomass with total neighbour biomass (writes `Table_S6_biomass_ancova.csv`, the file name of an earlier numbering) |
| `candidate_pool.R`, `table_s5_biomass_taxa.R` | qmd 3.2.2 | eligible taxa; TITAN2 indicator table (exploratory, pot level), the source of Table S7 (file names kept from an earlier numbering: the module writes `Table_S5_biomass_associated_taxa.csv`) |
| `trait_models_v3.R`, `species_traits.csv`, `weeds_wheat_VPM2_WP_S3.tre`, `table3_robustness.py` | qmd 3.3 | Table 2, Notes S3 |
| `main_figures.py` | notebook | Figs 1–4 (Times New Roman) |
| `si_figures.py`, `fig_alpha_diversity.py` | notebook | Figs S5–S7 |
| `da_consensus.py`, `da_core_check.py`, `bacteria_dotplot_revised.py`, `fungi_dotplot_cell.py` | notebook 3.1.2 | Fig. 2, Table S3, Methods S3 tests |
| `pooled_check.R`, `vs_alone.R`, `control_size.py` | notebook 3.1.2 | control-size and wheat-alone checks (Methods S3, Table S3) |
| `assembly_check.R` | stand-alone (`Rscript assembly_check.R <matrix_dir> <rpoB_tree.nwk> <out_dir>`) | turnover with betapart (Methods S2) and bacterial βNTI with iCAMP (Fig. S7) |
| `core_vs_indicators.py` | notebook 3.1.2 | overlap of Fig. 2 ASVs with TITAN2 indicators (Results 3.2.2) |
| `biomass_tests.py` | notebook 3.2.1 | log-scale ANOVA, between-soil tests (Fig. 3), legume test, Table S4 |
| `build_ga.py` | notebook | graphical abstract |
