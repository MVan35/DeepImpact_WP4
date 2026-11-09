## ===========================================================================
## Differential-abundance co-detection -- the published definition, in one place
##
## "Differentially enriched" means appearing in the UNION of three methods, each against the
## `wheat x wheat` control, each with its own rule:
##
##   DESeq2         |log2FC| > 2  AND unadjusted pvalue < 0.05
##   edgeR          |logFC|   > 2  AND unadjusted PValue < 0.01
##   metagenomeSeq  adjPvalues < 0.05
##
## ===========================================================================

## ---------------------------------------------------------------------------
## Resolving which copy of a DA export to read
##
## DA_SEARCH is the ordered list of roots to look under. 
## ---------------------------------------------------------------------------


## `control` defaults to CTRL, which da_codetection() sets with <<- before
da_resolve <- function(subdir, fname,
                       control = if (exists("CTRL", inherits = TRUE))
                                   get("CTRL", inherits = TRUE)
                                 else "wheat x wheat",
                       roots = NULL,
                       quiet = FALSE) {
  if (is.null(roots))
    roots <- if (exists("DA_SEARCH", inherits = TRUE) &&
                 !is.null(get("DA_SEARCH", inherits = TRUE)))
               get("DA_SEARCH", inherits = TRUE) else RES_DIR
  cand <- file.path(roots, subdir, fname)
  cand <- cand[file.exists(cand)]
  if (!length(cand)) {
    warning("no copy of ", fname, " under: ", paste(roots, collapse = ", "),
            call. = FALSE)
    return(NA_character_)
  }
  info <- lapply(cand, function(f) {
    d <- utils::read.csv(f, check.names = FALSE, nrows = 20000)
    has <- FALSE
    if ("comparison" %in% names(d)) {
      p <- do.call(rbind, strsplit(as.character(d$comparison), " vs ", fixed = TRUE))
      p <- gsub(".", " ", p, fixed = TRUE)
      has <- control %in% p[, 1] || control %in% p[, 2]
    } else has <- TRUE          # no comparison column -> nothing to check
    list(f = f, rows = nrow(d), has = has)
  })
  good <- Filter(function(x) x$has, info)
  if (!length(good)) {
    warning(sprintf(paste0("%s: no copy contains the control '%s'. Checked:\n  %s\n",
                           "The direction-aware filter will return an empty set ",
                           "without erroring, so this is reported rather than ",
                           "left silent."),
                    fname, control,
                    paste(sprintf("%s (%d rows)", sapply(info, `[[`, "f"),
                                  sapply(info, `[[`, "rows")), collapse = "\n  ")),
            call. = FALSE)
    return(NA_character_)
  }
  pick <- good[[which.max(sapply(good, `[[`, "rows"))]]
  if (!quiet && length(info) > 1)
    cat(sprintf("  %-52s <- %s (%d rows; %d candidate%s)\n", fname,
                basename(dirname(dirname(pick$f))), pick$rows, length(info),
                if (length(info) > 1) "s" else ""))
  pick$f
}

split_comparison <- function(d) {
  parts <- strsplit(as.character(d$comparison), " vs ", fixed = TRUE)
  d$species1 <- gsub(".", " ", vapply(parts, `[`, "", 1), fixed = TRUE)
  d$species2 <- gsub(".", " ", vapply(parts, `[`, "", 2), fixed = TRUE)
  d
}

## GUARD: if the control label is absent from a file, the direction-aware filter
## matches nothing and returns an EMPTY set without erroring. 
check_control <- function(d, label) {
  ok <- CTRL %in% d$species1 || CTRL %in% d$species2
  if (!ok)
    warning(sprintf("%s: control '%s' appears on NEITHER side of `comparison`. ",
                    label, CTRL),
            sprintf("File holds %d rows with species1 = {%s}. The filter will ",
                    nrow(d), paste(unique(d$species1), collapse = ", ")),
            "return an empty set. This is a wrong or truncated export, not a null result.",
            call. = FALSE)
  ok
}

da_deseq <- function(soil_tag) {
  f <- da_resolve("DESeq2_species", sprintf("DESeq_bacteria_%s_species_vs_wheat_x_wheat_pval.csv", soil_tag))
  if (is.na(f)) return(character(0))
  d <- read.csv(f, check.names = FALSE)
  if (!"taxa" %in% names(d)) {
    names(d) <- c("cluster", "baseMean", "log2FoldChange", "lfcSE", "stat",
                  "pvalue", "padj", "taxa", "species1", "species2")
  }
  d$cluster <- d$taxa
  d <- split_comparison(d)
  check_control(d, basename(f))
  keep <- (d$species1 == CTRL & d$log2FoldChange < -2 & d$pvalue < 0.05) |
          (d$species2 == CTRL & d$log2FoldChange >  2 & d$pvalue < 0.05)
  unique(as.character(d$cluster[which(keep)]))
}

da_edger <- function(soil_tag) {
  f <- da_resolve("edgeR_species", sprintf("edgeR_bacteria_%s_vs_wheat_x_wheat.csv", soil_tag))
  if (is.na(f)) return(character(0))
  d <- read.csv(f, check.names = FALSE)
  d$cluster <- if ("taxa" %in% names(d)) d$taxa else d[[1]]
  d <- split_comparison(d)
  check_control(d, basename(f))
  keep <- (d$species1 == CTRL & d$logFC < -2 & d$PValue < 0.01) |
          (d$species2 == CTRL & d$logFC >  2 & d$PValue < 0.01)
  unique(as.character(d$cluster[which(keep)]))
}

da_metaseq <- function(soil_tag) {
  f <- da_resolve("metagenomeSeq_species", sprintf("MetaSeq_bacteria_%s_CSS_wheat_x_wheat.csv", soil_tag))
  if (is.na(f)) return(character(0))
  d <- read.csv(f, check.names = FALSE)
  d$cluster <- if ("taxa" %in% names(d)) d$taxa else d[[1]]
  unique(as.character(d$cluster[which(d$adjPvalues < 0.05)]))
}

SOILTAG <- c(cereal = "Resistant", mixed = "Permissive")


da_codetection <- function(titan_dir, res_dir, control = "wheat x wheat",
                           cp_min = 0.35, purity_min = 0.95,
                           reliability_min = 0.95,
                           soiltag = c(cereal = "Resistant", mixed = "Permissive"),
                           soilkey = c(cereal = "soilA", mixed = "soilB")) {
  CTRL <<- control          # the loaders read this
  RES_DIR <<- res_dir
  counts <- list(); parts <- list()
  for (s in names(soiltag)) {
    tag <- soiltag[[s]]
    p <- list(DESeq2 = da_deseq(tag), edgeR = da_edger(tag),
              metagenomeSeq = da_metaseq(tag))
    parts[[s]] <- p
    ## the high-biomass (declining) threshold set for this soil -- the 15 and 25
    hb <- titan_flagged(file.path(titan_dir,
            sprintf("bacteria_titan_%s_results_v3.csv", soilkey[[s]])),
            cp_min = cp_min, purity_min = purity_min,
            reliability_min = reliability_min)[["z-"]]
    uni <- unique(unlist(p))
    ov <- intersect(hb, uni)
    counts[[s]] <- data.frame(
      soil = s, n_hb = length(hb),
      n_deseq = length(p$DESeq2), n_edger = length(p$edgeR),
      n_metaseq = length(p$metagenomeSeq), n_union = length(uni),
      overlap = length(ov), pct = 100 * length(ov) / length(hb),
      ov_deseq = length(intersect(hb, p$DESeq2)),
      ov_edger = length(intersect(hb, p$edgeR)),
      ov_metaseq = length(intersect(hb, p$metagenomeSeq)),
      taxa = paste(sort(ov), collapse = "; "),
      stringsAsFactors = FALSE)
  }
  ## unique contribution: taxa this method finds and no other. 
  uniq <- do.call(rbind, lapply(names(parts), function(s) {
    hb <- titan_flagged(file.path(titan_dir,
            sprintf("bacteria_titan_%s_results_v3.csv", soilkey[[s]])),
            cp_min = cp_min, purity_min = purity_min,
            reliability_min = reliability_min)[["z-"]]
    p <- parts[[s]]
    do.call(rbind, lapply(names(p), function(m) {
      others <- unique(unlist(p[setdiff(names(p), m)]))
      data.frame(soil = s, method = m,
                 unique_overlap = length(intersect(hb, setdiff(p[[m]], others))))
    }))
  }))
  list(counts = do.call(rbind, counts), unique = uniq,
       provenance = da_provenance(res_dir, control, soiltag), parts = parts)
}

da_provenance <- function(res_dir, control = "wheat x wheat",
                          soiltag = c(cereal = "Resistant", mixed = "Permissive")) {
  do.call(rbind, lapply(names(soiltag), function(s) {
    tag <- soiltag[[s]]
    do.call(rbind, lapply(c("DESeq2", "edgeR", "metagenomeSeq"), function(m) {
      ## Resolve exactly as the loaders do
      f <- switch(m,
        DESeq2 = da_resolve("DESeq2_species",
                 sprintf("DESeq_bacteria_%s_species_vs_wheat_x_wheat_pval.csv", tag),
                 quiet = TRUE),
        edgeR = da_resolve("edgeR_species",
                sprintf("edgeR_bacteria_%s_vs_wheat_x_wheat.csv", tag), quiet = TRUE),
        metagenomeSeq = da_resolve("metagenomeSeq_species",
                        sprintf("MetaSeq_bacteria_%s_CSS_wheat_x_wheat.csv", tag),
                        quiet = TRUE))
      if (is.na(f) || !file.exists(f))
        return(data.frame(soil = s, method = m, root = NA_character_,
                          rows = NA_integer_, control_side = "FILE MISSING",
                          stringsAsFactors = FALSE))
      d <- read.csv(f, check.names = FALSE)
      root <- basename(dirname(dirname(f)))
      if (!"comparison" %in% names(d))
        return(data.frame(soil = s, method = m, root = root, rows = nrow(d),
                          control_side = "n/a (no comparison column)",
                          stringsAsFactors = FALSE))
      d <- split_comparison(d)
      side <- if (control %in% d$species1) "species1"
              else if (control %in% d$species2) "species2"
              else "ABSENT -- BAD EXPORT"
      data.frame(soil = s, method = m, root = root, rows = nrow(d), control_side = side,
                 stringsAsFactors = FALSE)
    }))
  }))
}
