# Analysis code

Code for *Neighbouring weed species influence the root microbiota and growth of wheat* (Vanhove *et al.*, New Phytologist, submitted).

## What you need

* **Data**: the Zenodo data record (DOI in the article's Data availability statement). Unpack it as `data/` next to this folder (`../data`), or set `PAPER_DATA` to its `data/` folder. Results are written to `../results`, or `../results_from_asv` with the switch below (`PAPER_OUT` overrides).
* **R ≥ 4.2** with Quarto or RStudio. CRAN: `vegan`, `dplyr`, `ggplot2`, `patchwork`, `car`, `ape`, `knitr`. Bioconductor: `phyloseq`, `edgeR`, `DESeq2`, `metagenomeSeq`.
  Optional regenerations: `TITAN2`, `NST`, `picante`, `betapart`, `iCAMP` (`V.PhyloMaker2` only to rebuild the supplied weed tree).
* **Python ≥ 3.10** with Jupyter: `numpy`, `pandas`, `scipy`, `matplotlib`, `statsmodels`, `openpyxl`, `Pillow`. `Rscript` on the PATH (or `RSCRIPT` set) for the control-size check.

## Run order

1. `Paper_analyses.qmd` (R): every statistical analysis, organised by the Results sections. Writes to `results/` (`PAPER_OUT`). About 7 min.
2. `Paper_figures.ipynb` (Python): every figure, the Python-side tables and the display items (`results/display_items/`). Reads the data folder and the outputs of step 1. About 3 min.

Both documents are self-contained apart from the files below. They locate this folder themselves (`PAPER_CODE` overrides).

| File | Used by |
|---|---|
| `validity_checks.R`, `robustness_and_tables.R`, `variance_partition_plot.R`, `table_s7_titan2.R`, `da_codetection.R`, `edger_export.R` | qmd, sourced in `load-modules` |
| `pooled_check.R`, `vs_alone.R` | notebook 3.1.2 (control-size check, Methods S3, Table S3) |
| `species_traits.csv` | qmd 3.3 and notebook 3.3 (Table 2, Table S1 traits) |
| `weeds_wheat_VPM2_WP_S3.tre` | qmd 3.3 (weed phylogeny, V.PhyloMaker2 scenario S3) |

## Options

The defaults reproduce the article. Set these environment variables to `"true"` to regenerate intermediate files (written to `results/`, compared with the deposited ones):

| Variable | Regenerates | Time |
|---|---|---|
| `PAPER_RUN_DA` | DESeq2 and metagenomeSeq exports (runs automatically when missing) | c. 5 min |
| `PAPER_RUN_NST`, `PAPER_RUN_BNTI` | taxonomic and phylogenetic NST, βNTI (seeded) | c. 4 h (βNTI prints no progress for c. 1 h) |
| `PAPER_RUN_ASSEMBLY` | betapart turnover and iCAMP βNTI (Methods S2) | several hours (4 parallel workers) |
| `PAPER_RUN_TITAN` (`PAPER_TITAN_CELLS` for a subset) | TITAN2 exports (seeded) | c. 8.5 h |

**From the ASV tables only.** Set `RUN_FROM_ASV <- TRUE` at the top of the `paths` chunk of the qmd and `RUN_FROM_ASV = True` in the Setup cell of the notebook (or the environment variable `PAPER_FROM_ASV=true` for both): every intermediate result above is regenerated from the ASV tables and the two ASV trees and then used (c. 16 h on a 2025 laptop: TITAN2 c. 8.5 h, NST and βNTI c. 4 h, iCAMP several hours). The deposited exports are not needed in that case.

`PAPER_DA = "regenerated"`, `PAPER_TITAN = "regenerated"` and `PAPER_NST = "original"` switch the inputs used.
The deposited differential-abundance and TITAN2 exports are the ones of the article. Re-exports under other package
versions differ (the article used DESeq2 v1.48, edgeR v4.6, metagenomeSeq v1.52). Under DESeq2 1.50 and edgeR 4.8 the
bacterial ASVs of Fig. 2 are all recovered, with additional ones, and 9 of the 10 fungal ASVs.

## Display items

Figs 1–4, Tables 1–2, Figs S4–S7 and Tables S3–S7 are produced by the two documents (the notebook's last cell lists them).
Figs S1–S3, Tables S1–S2 and Methods S1 are not produced by this code.
