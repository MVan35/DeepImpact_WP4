## ===========================================================================
## Part A — CLR / Aitchison robustness  (helper section 3.3)
## Part B — PERMANOVA table built from the adonis2 objects, not hardcoded
##
## Drop next to Analyses_bacteria_v1.2.Rmd and source(). Entry point is a
## phyloseq object, so it uses physeq_bacteria / physeq_fungi directly.
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


## ═══════════════════════════════════════════════════════════════════════
## PART A — CLR / Aitchison robustness
## ═══════════════════════════════════════════════════════════════════════
##
## THE OBJECTION a referee will raise:
##   "You rarefied. McMurdie & Holmes (2014) call that inadmissible, and
##    Bray-Curtis on proportions ignores the compositional nature of the
##    data (Gloor et al. 2017)."
##
## THE ANSWER: McMurdie & Holmes is about DIFFERENTIAL ABUNDANCE testing of
## individual taxa, where discarding reads costs power. For beta-diversity and
## PERMANOVA, equalising depth is defensible (Weiss et al. 2017). But arguing
## is weaker than showing: if the variance partition is stable across
## normalisations, the objection is disposed of in one sentence.
##
## WHAT THIS CAN AND CANNOT TEST, given you no longer have the raw objects:
##   * CAN test  — the compositional GEOMETRY: does replacing Bray-Curtis on
##     proportions with Aitchison distance (Euclidean on CLR) change the
##     partition? This is the Gloor et al. objection, and it is the one that
##     actually bears on your conclusions.
##   * CANNOT test — the rarefying DECISION itself. That needs the
##     pre-rarefaction counts. Say so explicitly rather than implying you
##     tested it; see the suggested wording at the bottom.
##
## Note on CLR from rarefied counts: rarefaction is itself a compositional
## operation (it fixes the library sum), so CLR of a rarefied table is a
## legitimate compositional representation. It is not equivalent to CLR of the
## raw counts, which is why the caveat above matters.
## ---------------------------------------------------------------------------

## zero replacement --------------------------------------------------------
## CLR needs log(x), so zeros must go. Two options:
##   "half_min"  x + min(nonzero)/2. Simple, deterministic, reproducible.
##   "cmultRepl" Bayesian-multiplicative replacement (zCompositions). More
##               principled; preserves ratios among non-zero parts. Use it if
##               the package is installed, since a referee may ask.
clr_transform <- function(M, zero_method = c("cmultRepl", "half_min")) {
  zero_method <- match.arg(zero_method)
  M <- as.matrix(M); storage.mode(M) <- "double"

  if (zero_method == "cmultRepl" && requireNamespace("zCompositions", quietly = TRUE)) {
    # cmultRepl expects samples x parts and returns closed (proportional) data
    Mp <- as.matrix(zCompositions::cmultRepl(M, label = 0, method = "CZM",
                                             output = "p-counts", z.warning = 1))
    used <- "cmultRepl (CZM)"
  } else {
    if (zero_method == "cmultRepl")
      message("zCompositions not installed; falling back to half-minimum replacement")
    Mp <- M + min(M[M > 0]) / 2
    used <- sprintf("half-minimum (+%.3g)", min(M[M > 0]) / 2)
  }
  clr <- t(apply(Mp, 1, function(z) log(z) - mean(log(z))))
  dimnames(clr) <- dimnames(M)
  attr(clr, "zero_replacement") <- used
  clr
}


## the comparison ----------------------------------------------------------
## ps        : phyloseq object as analysed (rarefied)
## label     : name for the output rows
## returns   : list(long = full adonis2 output for each normalisation,
##                  wide = the three main terms side by side,
##                  max_shift = largest change in R2 percentage points)
compare_normalisations <- function(ps, label,
                                   soil_var = "soil", species_var = "species",
                                   permutations = 999, seed = 1,
                                   zero_method = "cmultRepl") {

  M <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) M <- t(M)
  M <- matrix(as.numeric(M), nrow(M), ncol(M), dimnames = dimnames(M))
  meta <- as(sample_data(ps), "data.frame")
  meta[[soil_var]]    <- factor(meta[[soil_var]])
  meta[[species_var]] <- factor(meta[[species_var]])

  f <- as.formula(sprintf("D ~ %s + %s + %s:%s",
                          soil_var, species_var, soil_var, species_var))

  run <- function(D, nm) {
    set.seed(seed)
    a <- adonis2(f, data = meta, permutations = permutations, by = "terms")
    d <- as.data.frame(a); d$term <- rownames(d)
    d$normalisation <- nm; d$R2_pct <- d$R2 * 100
    d[, c("normalisation", "term", "Df", "R2", "R2_pct", "F", "Pr(>F)")]
  }

  ## (1) as published: Bray-Curtis on the rarefied counts
  D <- vegdist(M, method = "bray")
  r1 <- run(D, "rarefied_bray")

  ## (2) Bray-Curtis on relative abundances. On an already-rarefied table this
  ##     is near-identical to (1) — included because it is the check people
  ##     name, and showing it is a no-op is itself informative.
  D <- vegdist(M / rowSums(M), method = "bray")
  r2 <- run(D, "relabund_bray")

  ## (3) Aitchison = Euclidean distance on CLR-transformed abundances
  clr <- clr_transform(M, zero_method = zero_method)
  zr  <- attr(clr, "zero_replacement")
  D <- dist(clr, method = "euclidean")
  r3 <- run(D, "clr_aitchison")

  long <- rbind(r1, r2, r3)
  long$kingdom <- label
  rownames(long) <- NULL

  keep <- c(soil_var, species_var, paste0(soil_var, ":", species_var), "Residual")
  sub  <- long[long$term %in% keep, c("term", "normalisation", "R2_pct")]
  wide <- reshape(sub, idvar = "term", timevar = "normalisation",
                  direction = "wide")
  names(wide) <- sub("R2_pct\\.", "", names(wide))
  wide$shift_pp <- with(wide, clr_aitchison - rarefied_bray)   # percentage points

  list(long = long, wide = wide,
       max_shift = max(abs(wide$shift_pp[wide$term != "Residual"])),
       zero_replacement = zr, label = label)
}


## how to report it --------------------------------------------------------
report_normalisation <- function(res) {
  cat("\n=====", res$label, "— normalisation sensitivity =====\n")
  print(res$wide, row.names = FALSE, digits = 4)
  cat(sprintf("\nzero replacement: %s\n", res$zero_replacement))
  cat(sprintf("largest shift in any main term: %.2f percentage points\n", res$max_shift))
  if (res$max_shift < 2) {
    cat("-> STABLE. Suggested supplementary sentence:\n")
    cat("   \"The variance partition was insensitive to the choice of\n",
        "   dissimilarity: repeating the PERMANOVA on centred log-ratio\n",
        "   transformed abundances with Aitchison distance changed each term's\n",
        sprintf("   R\u00b2 by less than %.1f percentage points (Table Sx).\"\n", ceiling(res$max_shift * 10) / 10)
    )
  } else {
    cat("-> NOT stable. Report both partitions and say which you interpret, and why.\n")
  }
  cat("\nCaveat that must accompany this: the pre-rarefaction count tables are\n",
      "no longer available, so this tests the compositional geometry (Gloor et\n",
      "al. 2017) but not the rarefying decision itself (McMurdie & Holmes 2014).\n")
  invisible(res)
}


## ═══════════════════════════════════════════════════════════════════════
## PART B — PERMANOVA table from the fitted objects
## ═══════════════════════════════════════════════════════════════════════
##
## Replaces the hardcoded blocks in the Rmd:
##     bacteria_data <- data.frame(Factor = ..., R2 = c(0.02623, ...), ...)
## Those literals are correct today and will silently go stale the first time
## the OTU tables, filtering, or rarefaction depth change. Reading them from
## the adonis2 object makes that impossible.
## ---------------------------------------------------------------------------

## fits : named list of adonis2 objects, e.g. list(Bacteria = bact_seq, Fungi = fung_seq)
## Pass the SEQUENTIAL fits (by = "terms") — marginal fits omit the main
## effects, so a partition cannot be built from them.
permanova_table <- function(fits,
                            term_labels = c(soil = "Soil of origin",
                                            species = "Neighbour identity",
                                            "soil:species" = "Soil \u00d7 neighbour",
                                            Residual = "Residual"),
                            digits = 3) {

  one <- function(fit, kingdom) {
    d <- as.data.frame(fit)
    d$term <- rownames(d)
    d <- d[d$term != "Total", ]

    # vegan names the interaction after the formula, and has used both
    # "Residual" and "Residuals" across versions — normalise both.
    d$term <- sub("^species:soil$", "soil:species", d$term)
    d$term <- sub("^Residuals$", "Residual", d$term)

    miss <- setdiff(names(term_labels), d$term)
    if (length(miss))
      stop("adonis2 output for ", kingdom, " has no row(s): ",
           paste(miss, collapse = ", "),
           ".\n  Rows present: ", paste(d$term, collapse = ", "),
           ".\n  If you passed a by = \"margin\" fit, refit with by = \"terms\".",
           call. = FALSE)

    d <- d[match(names(term_labels), d$term), ]
    data.frame(
      Kingdom  = kingdom,
      Term     = unname(term_labels),
      Df       = d$Df,
      SumOfSqs = round(d$SumOfSqs, digits),
      R2       = round(d$R2, 5),
      R2_pct   = round(d$R2 * 100, 2),
      F        = round(d$F, 3),
      p        = d$`Pr(>F)`,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, Map(one, fits, names(fits)))
  rownames(out) <- NULL

  # partition must sum to 100% — a sanity check the hardcoded version lacked
  for (kg in unique(out$Kingdom)) {
    tot <- sum(out$R2_pct[out$Kingdom == kg])
    if (abs(tot - 100) > 0.5)
      warning(kg, ": R2 sums to ", round(tot, 2), "%, not 100%. ",
              "Did you pass a marginal fit, or omit the Residual row?",
              call. = FALSE)
  }
  out
}


## formatted for the manuscript -------------------------------------------
## p-values as "< 0.001" rather than 0.001, which is what adonis2 returns at
## 999 permutations and is not the same claim.
permanova_table_pretty <- function(tab, nperm = 999, dash = "\u2014") {
  floor_p <- 1 / (nperm + 1)
  # The Residual row has no F and no p — print an em dash, not NA.
  tab$p_fmt  <- ifelse(is.na(tab$p), dash,
                ifelse(tab$p <= floor_p, sprintf("< %.3f", floor_p),
                                          sprintf("%.3f", tab$p)))
  tab$F_fmt  <- ifelse(is.na(tab$F), dash, sprintf("%.2f", tab$F))
  tab$R2_fmt <- sprintf("%.1f", tab$R2_pct)
  tab[, c("Kingdom", "Term", "Df", "R2_fmt", "F_fmt", "p_fmt")]
}


## effect-size categories, if you want to keep that column ------------------
## Kept from the original Rmd, but note in the caption that these thresholds
## are conventional and not calibrated for dissimilarity-based R2 — a
## Bray-Curtis R2 is not comparable to a linear-model R2.
categorize_effect <- function(r2) {
  cut(r2, breaks = c(-Inf, 0.01, 0.06, 0.14, Inf),
      labels = c("Negligible", "Small", "Medium", "Large"))
}
