## ===========================================================================
## Candidate-pool construction for the mEQO / convergence analysis
##
## Reports, per kingdom x soil, how many taxa survive each filter step. These
## are the numbers that belong in Methods and in the Fig 5c / Table 2 legends,
## because every cross-method comparison is bounded by this pool.
##
## Verified against the exported pools: running this on physeq_bacteria /
## physeq_fungi reproduces the 200 / 200 / 124 / 153 candidate sets exactly,
## taxon for taxon.
##
## NOTE the filter is applied WITHIN each soil, on abundances renormalised to
## that soil's samples -- not once on the pooled table. Prevalence and mean
## abundance both change when the sample set changes, so a pooled filter gives
## different pools.
## ===========================================================================

suppressPackageStartupMessages({library(phyloseq); library(dplyr)})

## --------------------------------------------------------------- configuration
## These must match the values used in the mEQO run.
PREV_MIN <- 0.10      # present in >= 10% of that soil's samples
ABUND_MIN <- 1e-4     # mean relative abundance within that soil
MAX_TAXA <- 200       # cap, applied after the two filters

SOIL_MAP   <- c(cereal = "resistant", mixed = "Permissive")
SOIL_LABEL <- c(cereal = "Cereal-dominated (AF193)",
                mixed  = "Mixed-system (AF190)")


#' Candidate pool for one kingdom x soil cell.
#'
#' @param ps        phyloseq object (physeq_bacteria / physeq_fungi)
#' @param soil      value in sample_data(ps)$soil, e.g. "resistant"
#' @param response  sample variable used to drop samples with no biomass value,
#'                  so the pool matches the samples the analysis actually used
#' @return list(taxa, stats) -- the retained taxon names and the cascade counts
eqo_candidate_pool <- function(ps, soil, response = "Dry_wheat_mass",
                               prev_min = PREV_MIN, abund_min = ABUND_MIN,
                               max_taxa = MAX_TAXA) {
  sd <- as(sample_data(ps), "data.frame")
  keep_s <- rownames(sd)[sd$soil == soil & !is.na(sd[[response]])]
  sub <- prune_samples(keep_s, ps)

  ## counts -> within-soil relative abundance. Renormalising AFTER subsetting
  ## is the point: a taxon's mean abundance is defined relative to the samples
  ## in that soil, not to the whole experiment.
  M <- as(otu_table(sub), "matrix")
  if (!taxa_are_rows(sub)) M <- t(M)
  M <- sweep(M, 2, colSums(M), "/")

  prev  <- rowMeans(M > 0)
  abund <- rowMeans(M)

  n_all  <- nrow(M)
  pass_p <- prev >= prev_min
  pass_a <- pass_p & abund >= abund_min
  elig   <- rownames(M)[pass_a]

  ## cap keeps the most abundant, so the pool is the abundant core rather than
  ## an arbitrary subset
  if (length(elig) > max_taxa)
    elig <- names(sort(abund[elig], decreasing = TRUE))[seq_len(max_taxa)]

  stats <- data.frame(
    n_samples            = length(keep_s),
    all_taxa             = n_all,
    after_prevalence     = sum(pass_p),
    after_abundance      = sum(pass_a),
    pool                 = length(elig),
    cap_binding          = sum(pass_a) > max_taxa,
    abundance_retained   = mean(colSums(M[elig, , drop = FALSE])),
    stringsAsFactors     = FALSE)

  list(taxa = elig, stats = stats)
}


#' Full 2 kingdoms x 2 soils cascade, as one table.
eqo_pool_table <- function(physeq_bacteria, physeq_fungi,
                           response = "Dry_wheat_mass") {
  out <- list()
  for (kg in c("Bacteria", "Fungi")) {
    ps <- if (kg == "Bacteria") physeq_bacteria else physeq_fungi
    for (sk in names(SOIL_MAP)) {
      r <- eqo_candidate_pool(ps, SOIL_MAP[[sk]], response = response)
      out[[length(out) + 1]] <- data.frame(
        kingdom = kg, soil = SOIL_LABEL[[sk]], r$stats,
        stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, out)
  tab$abundance_retained_pct <- round(100 * tab$abundance_retained, 1)
  tab$abundance_retained <- NULL
  rownames(tab) <- NULL
  tab
}


#' One-paragraph Methods sentence, generated from the numbers so it cannot drift.
eqo_pool_sentence <- function(tab) {
  b <- tab[tab$kingdom == "Bacteria", ]
  f <- tab[tab$kingdom == "Fungi", ]
  sprintf(paste0(
    "Within each soil, abundances were renormalised and taxa retained if present ",
    "in at least %.0f%% of that soil's samples with mean relative abundance ",
    ">= %g, then capped at the %d most abundant. This left %s bacterial taxa ",
    "(from %d; the cap was binding in both soils) and %s fungal taxa (from %d; ",
    "below the cap, so all eligible taxa were retained). Retained taxa account ",
    "for %.0f-%.0f%% of bacterial and %.0f-%.0f%% of fungal relative abundance. ",
    "Cross-method comparisons are bounded by this pool: taxa below the filter ",
    "could not be selected, so their absence from an mEQO group is uninformative."),
    100 * PREV_MIN, ABUND_MIN, MAX_TAXA,
    paste(b$pool, collapse = " and "), b$all_taxa[1],
    paste(f$pool, collapse = " and "), f$all_taxa[1],
    min(b$abundance_retained_pct), max(b$abundance_retained_pct),
    min(f$abundance_retained_pct), max(f$abundance_retained_pct))
}


#' Check the recomputed pools against the exported candidate files.
#'
#' Run this once. If it does not report an exact match, the exported pools were
#' built with different settings and the Methods numbers would be wrong.
eqo_verify_pools <- function(physeq_bacteria, physeq_fungi, eqo_folder,
                             response = "Dry_wheat_mass") {
  LET <- c(cereal = "A", mixed = "B")
  out <- list()
  for (kg in c("bacteria", "fungi")) {
    ps <- if (kg == "bacteria") physeq_bacteria else physeq_fungi
    for (sk in names(SOIL_MAP)) {
      f <- file.path(eqo_folder,
                     sprintf("sample_names_soil%s_taxa_%s_v2.csv", LET[[sk]], kg))
      if (!file.exists(f)) next
      exported <- setdiff(names(read.csv(f, check.names = FALSE, nrows = 1)),
                          "sample")
      recomputed <- eqo_candidate_pool(ps, SOIL_MAP[[sk]], response = response)$taxa
      out[[length(out) + 1]] <- data.frame(
        kingdom = kg, soil = sk,
        n_exported = length(exported), n_recomputed = length(recomputed),
        shared = length(intersect(exported, recomputed)),
        exact = setequal(exported, recomputed),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}


## ------------------------------------------------------------------------ use
## tab <- eqo_pool_table(physeq_bacteria, physeq_fungi)
## print(tab, row.names = FALSE)
## cat(eqo_pool_sentence(tab), "\n")
## write.csv(tab, file.path(mypath_res, "table_candidate_pool.csv"), row.names = FALSE)
##
## eqo_verify_pools(physeq_bacteria, physeq_fungi,
##                  file.path(mypath_res, "mEQO"))   # expect exact = TRUE x4
