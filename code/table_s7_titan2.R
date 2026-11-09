## ===========================================================================
## table_s7_titan2.R -- Table S7: TITAN2 indicators of wheat biomass
## (exploratory, pot level; Methods 2.8, Methods S3)
##
## Replaces table_s5_biomass_taxa.R, which also carried the mEQO columns
## (mEQO selection, 200-taxon candidate-pool cap, r_alone). One row per taxon
## flagged by TITAN2 (purity and reliability >= 0.95, change point >= CP_MIN):
##
##   kingdom, soil, ASV, Genus       identity (Genus from the deposited taxonomy)
##   TITAN2_direction                z+ (increases above the change point) / z-
##   change_point_g (IndVal-based), z_change_point_g (z-score-based), purity, reliability: from the TITAN2 export
##   in_eligible_set                 passes prevalence >= PREV_MIN and mean
##                                   relative abundance >= ABUND_MIN in that soil
##   spearman_rho_biomass            Spearman rho with wheat biomass, pot level
##
## Usage:
##   source("table_s7_titan2.R")
##   s7 <- table_s7(physeq_bacteria, physeq_fungi, titan_dir = TITAN_DIR, matrix_dir = MATRIX_DIR)
## ===========================================================================
suppressPackageStartupMessages(library(phyloseq))

## TITAN2 taxa per direction (used by the qmd chunk `titan2`)
titan_flagged <- function(titan_csv, cp_min = 0.35, purity_min = 0.95, reliability_min = 0.95) {
  d <- read.csv(titan_csv, check.names = FALSE)
  if (!("taxa" %in% names(d)) && !nzchar(names(d)[1])) names(d)[1] <- "taxa"   # unnamed first column in the exports
  id <- if ("taxa" %in% names(d)) "taxa" else names(d)[1]
  d <- d[d$purity >= purity_min & d$reliability >= reliability_min & d$`ienv.cp` >= cp_min, ]
  list(`z-` = as.character(d[[id]][d$maxgrp == 1]), `z+` = as.character(d[[id]][d$maxgrp == 2]))
}

table_s7_cell <- function(ps, soil, soil_label, kingdom, titan_csv, genus,
                          cp_min, purity_min, reliability_min, prev_min, abund_min,
                          response = "Dry_wheat_mass") {
  sub <- prune_samples(sample_data(ps)$soil == soil, ps)
  M <- as(otu_table(sub), "matrix"); if (taxa_are_rows(sub)) M <- t(M)
  M <- M / rowSums(M)
  y <- as.numeric(sample_data(sub)[[response]]); keep <- !is.na(y); M <- M[keep, , drop = FALSE]; y <- y[keep]
  eligible <- colnames(M)[colMeans(M > 0) >= prev_min & colMeans(M) >= abund_min]

  d <- read.csv(titan_csv, check.names = FALSE)
  if (!("taxa" %in% names(d)) && !nzchar(names(d)[1])) names(d)[1] <- "taxa"
  id <- if ("taxa" %in% names(d)) "taxa" else names(d)[1]
  d <- d[d$purity >= purity_min & d$reliability >= reliability_min & d$`ienv.cp` >= cp_min, ]
  if (!nrow(d)) return(NULL)
  tx <- as.character(d[[id]])
  rho <- vapply(tx, function(t) if (t %in% colnames(M)) suppressWarnings(cor(M[, t], y, method = "spearman")) else NA_real_, numeric(1))
  out <- data.frame(kingdom = kingdom, soil = soil_label, ASV = tx, Genus = unname(genus[tx]),
                    TITAN2_direction = ifelse(d$maxgrp == 2, "z+", "z-"),
                    change_point_g = d$`ienv.cp`, purity = d$purity, reliability = d$reliability,
                    z_change_point_g = if ("zenv.cp" %in% names(d)) d$`zenv.cp` else NA_real_,
                    in_eligible_set = ifelse(tx %in% eligible, "yes", "no"),
                    spearman_rho_biomass = unname(rho), stringsAsFactors = FALSE, row.names = NULL)
  out[order(out$TITAN2_direction, out$ASV), ]
}

table_s7 <- function(physeq_bacteria, physeq_fungi, titan_dir, matrix_dir,
                     soils  = c(cereal = "resistant", mixed = "Permissive"),
                     labels = c(cereal = "Cereal-dominated (AF193)", mixed = "Mixed-system (AF190)"),
                     keys   = c(cereal = "soilA", mixed = "soilB"),
                     cp_min = 0.35, purity_min = 0.95, reliability_min = 0.95,
                     prev_min = 0.10, abund_min = 1e-4) {
  ps_of  <- list(bacteria = physeq_bacteria, fungi = physeq_fungi)
  tax_of <- c(bacteria = "tax_bacteria_R_1314ASVs.csv", fungi = "tax_fungi_R_647ASVs.csv")
  out <- list()
  for (kg in names(ps_of)) {
    tx <- read.delim(file.path(matrix_dir, tax_of[[kg]]), row.names = 1, check.names = FALSE)
    genus <- setNames(as.character(tx$Genus), rownames(tx))
    for (sl in names(soils))
      out[[length(out) + 1]] <- table_s7_cell(ps_of[[kg]], soils[[sl]], labels[[sl]], kg,
        file.path(titan_dir, sprintf("%s_titan_%s_results_v3.csv", kg, keys[[sl]])), genus,
        cp_min, purity_min, reliability_min, prev_min, abund_min)
  }
  s7 <- do.call(rbind, out); rownames(s7) <- NULL
  s7
}
