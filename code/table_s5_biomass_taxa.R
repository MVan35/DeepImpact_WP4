## ===========================================================================
## Table S5 -- Biomass-associated taxa by method
##
## REPLACES the old "multi-method microbial indicators" table, which was built
## around a >=2-method co-detection rule. That rule is retired: it is not a
## usable confidence criterion when one method returns 10 taxa and the other
## returns 7-61 from pools of 124-200 (see Table S10 for the chance
## expectations).
##
## The replacement is a per-taxon audit. One row per taxon flagged by EITHER
## method, with four columns that let a reader check any claim in section 3.2.2:
##
##   TITAN2_direction  z+ / z- / "-"     which change-point direction, if any
##   mEQO              selected / "-"    whether it is in the Nmax = 10 group
##   in_candidate_pool_capped            could mEQO have selected it? (200 cap)
##   in_eligible_set_uncapped            does it pass the two filters at all?
##
## The two pool columns matter and are NOT redundant. A TITAN2 indicator with
## capped = "no", uncapped = "yes" was excluded by the 200-taxon truncation, not
## by biology -- mEQO could not have selected it, so its absence from the mEQO
## column is uninformative. A taxon with uncapped = "no" fails the prevalence or
## abundance filter and is outside the comparison entirely. Blank cells in the
## old table conflated these two situations with genuine disagreement.
##
## Usage:
##   source("table_s5_biomass_taxa.R")
##   s5 <- table_s5(physeq_bacteria, physeq_fungi,
##                  titan_dir = "paper_1/data/results_v4/TITAN_Full_Dataset",
##                  meqo_dir  = "paper_1/data/results_v4/mEQO")
##   write.csv(s5, "Table_S5_biomass_associated_taxa.csv", row.names = FALSE)
##
## NOTE on the TITAN2 directory: paper_1/data/ and paper_1/data/results_v4/ both
## contain files named *_titan_soil[AB]_results_v3.csv and they are DIFFERENT.
## results_v4/ is the one that reproduces the flagged counts in the manuscript
## (24/15 bacterial z+/z- in cereal-dominated, 36/25 in mixed-system). Point
## titan_dir at results_v4/.
## ===========================================================================

library(phyloseq)

PREV_MIN <- 0.10; ABUND_MIN <- 1e-4; MAX_TAXA <- 200
CP_MIN   <- 0.35; PURITY_MIN <- 0.95; RELIABILITY_MIN <- 0.95

## ---------------------------------------------------------------------------
## 1. Abundance matrix, and the two pool definitions
## ---------------------------------------------------------------------------
soil_matrix <- function(ps, soil) {
  sub <- prune_samples(sample_data(ps)$soil == soil, ps)
  M <- as(otu_table(sub), "matrix"); if (taxa_are_rows(sub)) M <- t(M)
  M <- M / rowSums(M)
  y <- as.numeric(sample_data(sub)$Dry_wheat_mass)
  keep <- !is.na(y)
  list(M = M[keep, , drop = FALSE], y = y[keep])
}

apply_filters <- function(M, cap = NULL) {
  M <- M[, colMeans(M > 0) >= PREV_MIN & colMeans(M) >= ABUND_MIN, drop = FALSE]
  if (!is.null(cap) && ncol(M) > cap)
    M <- M[, order(colMeans(M), decreasing = TRUE)[seq_len(cap)], drop = FALSE]
  colnames(M)
}

## ---------------------------------------------------------------------------
## 2. TITAN2 indicators, split by direction
## ---------------------------------------------------------------------------
## TITAN2 taxa per direction with explicit thresholds (used by the qmd chunk `titan2` and by da_codetection.R)
titan_flagged <- function(titan_csv, cp_min = 0.35, purity_min = 0.95, reliability_min = 0.95) {
  d <- read.csv(titan_csv, check.names = FALSE)
  if (!("taxa" %in% names(d)) && !nzchar(names(d)[1])) names(d)[1] <- "taxa"
  id <- if ("taxa" %in% names(d)) "taxa" else names(d)[1]
  d <- d[d$purity >= purity_min & d$reliability >= reliability_min & d$`ienv.cp` >= cp_min, ]
  list(`z-` = as.character(d[[id]][d$maxgrp == 1]), `z+` = as.character(d[[id]][d$maxgrp == 2]))
}

titan_by_direction <- function(titan_csv) {
  d <- read.csv(titan_csv, check.names = FALSE)
  # taxon names sit in an UNNAMED first column in these exports; reading them
  # without this line silently returns an empty set rather than erroring
  if (!("taxa" %in% names(d)) && !nzchar(names(d)[1])) names(d)[1] <- "taxa"
  id <- if ("taxa" %in% names(d)) "taxa" else names(d)[1]
  d <- d[d$purity >= PURITY_MIN & d$reliability >= RELIABILITY_MIN &
         d$`ienv.cp` >= CP_MIN, ]
  list(`z+` = as.character(d[[id]][d$maxgrp == 2]),    # increases above cp
       `z-` = as.character(d[[id]][d$maxgrp == 1]))    # declines above cp
}

## ---------------------------------------------------------------------------
## 3. One cell of the table
## ---------------------------------------------------------------------------
table_s5_cell <- function(ps, soil, soil_label, kingdom, titan_csv, meqo_csv) {
  sm <- soil_matrix(ps, soil)
  pool_uncapped <- apply_filters(sm$M, cap = NULL)
  pool_capped   <- apply_filters(sm$M, cap = MAX_TAXA)

  # Spearman rho is computed on the UNCAPPED eligible set, so every TITAN2
  # indicator that passes the filters gets a value even if the cap excluded it
  rho <- sapply(pool_uncapped,
                function(tx) suppressWarnings(
                  cor(sm$M[, tx], sm$y, method = "spearman")))

  tf   <- titan_by_direction(titan_csv)
  ## the mEQO export is optional (mEQO is not part of the paper); without it the table holds TITAN2 taxa only
  meqo <- if (!is.null(meqo_csv) && file.exists(meqo_csv)) read.csv(meqo_csv) else
            data.frame(taxon = character(0), Genus = character(0), r_alone = numeric(0))
  sel  <- as.character(meqo$taxon)
  gen  <- setNames(as.character(meqo$Genus), sel)
  ra   <- setNames(meqo$r_alone, sel)

  taxa <- union(unlist(tf, use.names = FALSE), sel)
  dir_of <- function(tx) if (tx %in% tf[["z+"]]) "z+" else
                         if (tx %in% tf[["z-"]]) "z-" else "-"

  out <- data.frame(
    kingdom = kingdom, soil = soil_label, ASV = taxa,
    Genus = ifelse(taxa %in% sel, gen[taxa], ""),
    TITAN2_direction = vapply(taxa, dir_of, character(1)),
    mEQO = ifelse(taxa %in% sel, "selected", "-"),
    in_candidate_pool_capped   = ifelse(taxa %in% pool_capped,   "yes", "no"),
    in_eligible_set_uncapped   = ifelse(taxa %in% pool_uncapped, "yes", "no"),
    spearman_rho_biomass = unname(rho[taxa]),
    r_alone_if_selected  = ifelse(taxa %in% sel, ra[taxa], NA_real_),
    stringsAsFactors = FALSE, row.names = NULL)

  # top decile of |rho| within this cell -- the claim "the most strongly
  # biomass-associated taxa" in Results 3.2.2 reads off this column
  thr <- quantile(abs(rho), 0.90, na.rm = TRUE)
  out$top_decile <- ifelse(!is.na(out$spearman_rho_biomass) &
                           abs(out$spearman_rho_biomass) >= thr, "yes", "no")
  out[order(out$TITAN2_direction, out$ASV), ]
}

## ---------------------------------------------------------------------------
## 4. All four cells
## ---------------------------------------------------------------------------
table_s5 <- function(physeq_bacteria, physeq_fungi, titan_dir, meqo_dir,
                     soils = c(cereal = "resistant", mixed = "Permissive"),
                     labels = c(cereal = "Cereal-dominated (AF193)",
                                mixed  = "Mixed-system (AF190)")) {
  ps_of  <- list(bacteria = physeq_bacteria, fungi = physeq_fungi)
  key_of <- c(cereal = "soilA", mixed = "soilB")
  out <- list()
  for (kg in names(ps_of)) for (sl in names(soils)) {
    out[[length(out) + 1]] <- table_s5_cell(
      ps_of[[kg]], soils[[sl]], labels[[sl]], kg,
      file.path(titan_dir, sprintf("%s_titan_%s_results_v3.csv",
                                   kg, key_of[[sl]])),
      file.path(meqo_dir,  sprintf("eqo_selected_%s_taxa_%s_v2.csv",
                                   key_of[[sl]], kg)))
  }
  do.call(rbind, out)
}

## ---------------------------------------------------------------------------
## 5. The numbers the Results text quotes
## ---------------------------------------------------------------------------
table_s5_summary <- function(s5) {
  agg <- do.call(rbind, lapply(split(s5, list(s5$kingdom, s5$soil), drop = TRUE),
    function(g) data.frame(
      kingdom = g$kingdom[1], soil = g$soil[1], rows = nrow(g),
      titan_flagged = sum(g$TITAN2_direction != "-"),
      meqo_selected = sum(g$mEQO == "selected"),
      both = sum(g$TITAN2_direction != "-" & g$mEQO == "selected"),
      titan_lost_to_cap = sum(g$TITAN2_direction != "-" &
                              g$in_candidate_pool_capped == "no" &
                              g$in_eligible_set_uncapped == "yes"),
      top_decile = sum(g$top_decile == "yes"))))
  rownames(agg) <- NULL
  cat("\n-- Table S5 summary --\n"); print(agg, row.names = FALSE)
  cat("\n  'both' is the count for Results 3.2.2. 'titan_lost_to_cap' is the\n",
      " number of TITAN2 indicators the 200-taxon truncation removed from the\n",
      " search space -- these are NOT disagreements between the methods.\n")
  shared <- s5[s5$TITAN2_direction != "-" & s5$mEQO == "selected", ]
  cat("\n-- taxa flagged by both methods --\n")
  print(shared[, c("kingdom", "soil", "ASV", "Genus", "TITAN2_direction",
                   "spearman_rho_biomass")], row.names = FALSE, digits = 3)
  invisible(agg)
}
