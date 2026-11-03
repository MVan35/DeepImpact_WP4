## ===========================================================================
## Validity checks for the PERMANOVA variance partition — every assumption
##
## One function per check, then run_validity_checks() to do all of them and
## return a tidy list of data frames. Entry point is a phyloseq object, so it
## drops into Analyses_bacteria_v1.2.Rmd using physeq_bacteria / physeq_fungi.
##
## Each returned table is a data frame, so when numbers change you re-run the
## chunk and the report text reads off the tables rather than being retyped.
##
## THE FIVE CHECKS
##   1. reproduce   sequential partition (what the manuscript reports)
##   2. ss_type     Type I vs Type II, and Type I with term order swapped
##                  -> is the partition an artefact of entry order?
##   3. dispersion  PERMDISP, both centroid and spatial-median
##                  -> is a "location" effect actually a spread effect?
##   4. replication design audit: does the pot column imply pseudoreplication?
##   5. normalisation  rarefied Bray-Curtis vs unrarefied relative abundance
##                     vs CLR/Aitchison  -> the compositional objection (3.3)
## ===========================================================================

## NOTE on coercing sample_data(). as.data.frame() on a phyloseq
## sample_data object returns an object that still carries the S4 class, and
## terms.formula() then rejects it:
##   Error in terms.formula(formula, data = data) :
##     'data' argument is of the wrong type
## as(x, "data.frame") strips it properly. Every coercion in this file uses
## that form; do not "simplify" it back.

library(phyloseq)
library(vegan)
library(dplyr)


## ---------------------------------------------------------------------------
## helper: phyloseq -> (samples x taxa) matrix + metadata, aligned
## ---------------------------------------------------------------------------
vc_matrix <- function(ps) {
  M <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) M <- t(M)
  list(M = M, meta = as(sample_data(ps), "data.frame"))
}

vc_factors <- function(meta, soil_var = "soil", species_var = "species") {
  meta[[soil_var]]    <- factor(meta[[soil_var]])
  meta[[species_var]] <- factor(meta[[species_var]])
  meta
}

## long -> wide, base R only (no tidyr dependency, so these functions load in a
## bare R session). Row order follows first appearance, not alphabetical.
vc_widen <- function(df, row, col, val) {
  rl <- unique(df[[row]]); cl <- unique(df[[col]])
  out <- data.frame(rl, stringsAsFactors = FALSE); names(out) <- row
  for (cc in cl) {
    idx <- match(rl, df[[row]][df[[col]] == cc])
    out[[cc]] <- df[[val]][df[[col]] == cc][idx]
  }
  out
}


## ---------------------------------------------------------------------------
## CHECK 1 + 2 — partition, SS type, and term order
## ---------------------------------------------------------------------------
## Type I  (by = "terms")  attributes shared variance to whichever term enters
##                         FIRST. This is what the manuscript reports.
##
## Type II (by = "margin") tests each term against the full model. IMPORTANT:
##                         vegan reports ONLY the interaction here, and drops the
##                         main effects entirely — for either `a * b` or the
##                         explicit `a + b + a:b` form. This is deliberate, not a
##                         quirk: a main effect that is contained in a
##                         higher-order term has no well-defined marginal test.
##                         So do NOT expect (or quote) a "marginal R2 for soil";
##                         there isn't one. The marginal row you do get confirms
##                         the interaction is order-independent.
##
## The order swap is therefore the decisive check on the main effects: if soil's
## R2 barely moves when species enters first, the published partition is not an
## ordering artefact. That comparison, not Type II, is what answers "is the soil
## term real or just first in line?".
## ---------------------------------------------------------------------------
vc_ss_type <- function(D, meta, soil_var = "soil", species_var = "species",
                       permutations = 999, seed = 1) {

  f_seq  <- as.formula(sprintf("D ~ %s + %s + %s:%s",
                               soil_var, species_var, soil_var, species_var))
  f_swap <- as.formula(sprintf("D ~ %s + %s + %s:%s",
                               species_var, soil_var, soil_var, species_var))

  set.seed(seed)
  a_seq  <- adonis2(f_seq,  data = meta, permutations = permutations, by = "terms")
  set.seed(seed)
  a_swap <- adonis2(f_swap, data = meta, permutations = permutations, by = "terms")
  set.seed(seed)
  a_marg <- adonis2(f_seq,  data = meta, permutations = permutations, by = "margin")

  tidy <- function(a, lab) {
    d <- as.data.frame(a); d$term <- rownames(d); d$model <- lab
    # the swapped formula names the interaction "species:soil"; rename so the
    # three models line up on one row instead of producing a phantom 4th term
    d$term[d$term == paste0(species_var, ":", soil_var)] <-
      paste0(soil_var, ":", species_var)
    d[, c("model", "term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
  }
  out <- rbind(tidy(a_seq,  "TypeI_soil_first"),
               tidy(a_swap, "TypeI_species_first"),
               tidy(a_marg, "TypeII_marginal"))
  out$R2_pct <- out$R2 * 100
  rownames(out) <- NULL
  attr(out, "fits") <- list(sequential = a_seq, swapped = a_swap, marginal = a_marg)
  out
}


## ---------------------------------------------------------------------------
## CHECK 3 — dispersion homogeneity (PERMDISP)
## ---------------------------------------------------------------------------
## PERMANOVA's H0 is "same centroid AND same dispersion". Rejecting tells you
## one of them is false. PERMDISP is what lets you say it was the centroid.
##
## betadisper()'s default is type = "median" (spatial median). type = "centroid"
## is more sensitive and can disagree — report which one you used, because a
## referee recomputing with the default will get a different number otherwise.
## ---------------------------------------------------------------------------
vc_dispersion <- function(D, meta, vars = c("soil", "species"),
                          permutations = 999, seed = 3) {
  rows <- list()
  for (v in vars) for (ty in c("median", "centroid")) {
    set.seed(seed)
    bd <- betadisper(D, meta[[v]], type = ty)
    pt <- permutest(bd, permutations = permutations)
    rows[[paste(v, ty)]] <- data.frame(
      term = v, type = ty,
      F = pt$tab$F[1], p = pt$tab$`Pr(>F)`[1],
      # numeric magnitude matters as much as the p-value: a "significant"
      # dispersion difference of 0.598 vs 0.606 is not driving anything
      min_dist = min(tapply(bd$distances, meta[[v]], mean)),
      max_dist = max(tapply(bd$distances, meta[[v]], mean)),
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$spread_ratio <- out$max_dist / out$min_dist
  rownames(out) <- NULL
  out
}


## ---------------------------------------------------------------------------
## CHECK 4 — design audit: is the pot column pseudoreplication?
## ---------------------------------------------------------------------------
## A referee seeing ~8 samples sharing a `pot` value will assume 8
## pseudoreplicates and conclude the residual df is inflated eightfold.
##
## The test: if every pot level maps to exactly ONE species x soil combination
## AND the number of pot levels equals the number of design cells, then `pot`
## is a TREATMENT CODE (redundant with species x soil) and the samples sharing
## it are independent replicate pots. If instead several samples share a pot
## level WITHIN a cell that has more pot levels than cells, it is a real
## physical container and you need it as a random effect / permutation stratum.
## ---------------------------------------------------------------------------
vc_replication <- function(meta, pot_var = "pot", soil_var = "soil",
                           species_var = "species") {
  if (!pot_var %in% colnames(meta))
    return(data.frame(check = "pot column absent", verdict = "n/a"))

  per_pot <- meta %>%
    group_by(.data[[pot_var]]) %>%
    summarise(n_samples = n(),
              n_species = n_distinct(.data[[species_var]]),
              n_soil    = n_distinct(.data[[soil_var]]), .groups = "drop")

  n_pot_levels <- nrow(per_pot)
  n_cells      <- meta %>% distinct(.data[[soil_var]], .data[[species_var]]) %>% nrow()
  mixed        <- sum(per_pot$n_species > 1 | per_pot$n_soil > 1)
  # is species x soil -> pot one-to-one?
  cells_multi_pot <- meta %>%
    group_by(.data[[soil_var]], .data[[species_var]]) %>%
    summarise(np = n_distinct(.data[[pot_var]]), .groups = "drop") %>%
    filter(np > 1) %>% nrow()

  treatment_code <- (mixed == 0 && n_pot_levels == n_cells && cells_multi_pot == 0)

  data.frame(
    n_samples        = nrow(meta),
    n_pot_levels     = n_pot_levels,
    n_design_cells   = n_cells,
    pots_mixing_cells = mixed,
    cells_with_multiple_pots = cells_multi_pot,
    median_samples_per_pot = median(per_pot$n_samples),
    verdict = if (treatment_code)
      "pot is a TREATMENT CODE (= species x soil); samples are independent replicates; free permutation OK"
    else
      "pot may be a PHYSICAL container: check before trusting residual df; consider strata = pot",
    stringsAsFactors = FALSE)
}


## ---------------------------------------------------------------------------
## CHECK 5 — §3.3 rarefying and the compositional objection
## ---------------------------------------------------------------------------
## The objection: McMurdie & Holmes (2014) call rarefying "inadmissible".
## The answer: their argument is about DIFFERENTIAL ABUNDANCE testing of
## individual taxa, where discarding reads costs power. For beta-diversity and
## PERMANOVA, equalising depth is defensible and standard (Weiss et al. 2017,
## Microbiome 5:27) because Bray-Curtis responds to library-size differences.
##
## The way to close it permanently is not to argue but to show the partition is
## stable across normalisations. This runs three:
##   (a) rarefied Bray-Curtis        — as published
##   (b) unrarefied relative abundance, Bray-Curtis  — no reads discarded
##   (c) CLR + Euclidean (Aitchison) — the compositional treatment
##       (Gloor et al. 2017, Front Microbiol 8:2224)
##
## If soil/species/interaction R2 move only a little across (a)-(c), one
## supplementary sentence disposes of the objection.
##
## ps_raw : phyloseq with RAW counts (pre-rarefaction)
## depth  : rarefaction depth used in the manuscript (7000 bact / 2000 fungi)
## ---------------------------------------------------------------------------
vc_normalisation <- function(ps_raw, depth, soil_var = "soil",
                             species_var = "species", permutations = 999,
                             seed = 1, clr_pseudocount = "auto") {

  f <- as.formula(sprintf("D ~ %s + %s + %s:%s",
                          soil_var, species_var, soil_var, species_var))
  run <- function(D, meta, lab) {
    set.seed(seed)
    a <- adonis2(f, data = meta, permutations = permutations, by = "terms")
    d <- as.data.frame(a); d$term <- rownames(d); d$normalisation <- lab
    d$R2_pct <- d$R2 * 100
    d[, c("normalisation", "term", "Df", "R2", "R2_pct", "F", "Pr(>F)")]
  }

  res <- list()

  ## (a) rarefied Bray-Curtis — as published
  set.seed(seed)
  ps_rar <- rarefy_even_depth(ps_raw, sample.size = depth,
                              rngseed = seed, replace = FALSE, verbose = FALSE)
  x <- vc_matrix(ps_rar); m <- vc_factors(x$meta, soil_var, species_var)
  res$rarefied <- run(vegdist(x$M, "bray"), m,
                      sprintf("rarefied_%d_bray", depth))

  ## (b) unrarefied relative abundance — nothing discarded
  ps_rel <- transform_sample_counts(ps_raw, function(z) if (sum(z) > 0) z/sum(z) else z)
  x2 <- vc_matrix(ps_rel); m2 <- vc_factors(x2$meta, soil_var, species_var)
  res$relative <- run(vegdist(x2$M, "bray"), m2, "unrarefied_relabund_bray")

  ## (c) CLR + Euclidean = Aitchison
  x3 <- vc_matrix(ps_raw); m3 <- vc_factors(x3$meta, soil_var, species_var)
  M <- x3$M
  # CLR needs zeros handled. Default: half the smallest non-zero value, which
  # is a simple and reproducible choice. zCompositions::cmultRepl() is the more
  # principled Bayesian-multiplicative alternative if it is installed.
  pc <- if (identical(clr_pseudocount, "auto")) min(M[M > 0]) / 2 else clr_pseudocount
  Mp <- M + pc
  clr <- t(apply(Mp, 1, function(z) log(z) - mean(log(z))))
  res$clr <- run(dist(clr, method = "euclidean"), m3, "clr_aitchison")

  out <- do.call(rbind, res)
  rownames(out) <- NULL

  # wide comparison of the three main terms — this is the table for the SI
  keep_terms <- c(soil_var, species_var, paste0(soil_var, ":", species_var))
  sub <- out[out$term %in% keep_terms, c("term", "normalisation", "R2_pct")]
  wide <- vc_widen(sub, row = "term", col = "normalisation", val = "R2_pct")

  list(long = out, wide = wide,
       clr_pseudocount = pc, depth = depth)
}


## ---------------------------------------------------------------------------
## RUN EVERYTHING
## ---------------------------------------------------------------------------
## ps_rarefied : the object the manuscript analysed (already rarefied)
## ps_raw      : raw-count object, for the normalisation check. If NULL that
##               check is skipped (you cannot un-rarefy).
## ---------------------------------------------------------------------------
run_validity_checks <- function(ps_rarefied, label,
                                ps_raw = NULL, depth = NULL,
                                soil_var = "soil", species_var = "species",
                                pot_var = "pot", permutations = 999, seed = 1) {

  x <- vc_matrix(ps_rarefied)
  meta <- vc_factors(x$meta, soil_var, species_var)
  D <- vegdist(x$M, method = "bray")

  message("[", label, "] 1-2. partition, SS type, term order ...")
  ss <- vc_ss_type(D, meta, soil_var, species_var, permutations, seed)

  message("[", label, "] 3. dispersion (PERMDISP) ...")
  disp <- vc_dispersion(D, meta, c(soil_var, species_var), permutations, seed)

  message("[", label, "] 4. design / replication audit ...")
  rep_ <- vc_replication(meta, pot_var, soil_var, species_var)

  norm <- NULL
  if (!is.null(ps_raw) && !is.null(depth)) {
    message("[", label, "] 5. normalisation sensitivity (rarefied / relabund / CLR) ...")
    norm <- vc_normalisation(ps_raw, depth, soil_var, species_var, permutations, seed)
  } else {
    message("[", label, "] 5. skipped (supply ps_raw + depth for the compositional check)")
  }

  # headline partition, for the report text and the abstract figure
  partition <- ss %>%
    filter(model == "TypeI_soil_first") %>%
    select(term, Df, R2_pct, F, `Pr(>F)`)

  list(label = label, partition = partition, ss_type = ss,
       dispersion = disp, replication = rep_, normalisation = norm,
       n_samples = nrow(x$M), n_taxa = ncol(x$M))
}


## ---------------------------------------------------------------------------
## Compact console report
## ---------------------------------------------------------------------------
print_validity <- function(v) {
  cat("\n=====", v$label, "—", v$n_samples, "samples x", v$n_taxa, "taxa =====\n")
  cat("\n-- Partition (Type I, soil first; what the manuscript reports) --\n")
  print(v$partition, row.names = FALSE, digits = 4)

  cat("\n-- SS type / term order (R2 %) --\n")
  sub <- v$ss_type[v$ss_type$term != "Total", c("term", "model", "R2_pct")]
  print(vc_widen(sub, "term", "model", "R2_pct"), row.names = FALSE, digits = 4)
  cat("  -> compare the two Type I columns: if they agree, entry order is not\n",
      "    driving the partition. The marginal column shrinking is expected\n",
      "    when the interaction is strong; report it, do not hide it.\n")

  cat("\n-- Dispersion (PERMDISP) --\n")
  print(v$dispersion, row.names = FALSE, digits = 4)
  cat("  -> spread_ratio near 1 means the groups have near-identical spread,\n",
      "    so a significant PERMANOVA is a centroid shift even if p is small.\n")

  cat("\n-- Design / replication --\n")
  print(t(v$replication), quote = FALSE)

  if (!is.null(v$normalisation)) {
    cat("\n-- Normalisation sensitivity (R2 %, section 3.3) --\n")
    print(v$normalisation$wide, row.names = FALSE, digits = 4)
    cat("  -> stable rows across the three columns close the rarefying /\n",
        "    compositional objection. CLR pseudocount used:",
        signif(v$normalisation$clr_pseudocount, 3), "\n")
  }
  invisible(v)
}
