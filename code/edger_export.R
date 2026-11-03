## edgeR: every neighbour vs the wheat x wheat control, one soil, one kingdom.
## Same model as Analyses_bacteria_v1.2.Rmd §EdgeR (L8055-8127): taxa present in >= 5
## samples of that soil, TMM, ~0 + species, estimateDisp, glmQLFit, glmQLFTest on
## "<neighbour> - wheat.x.wheat", rows with raw PValue < 0.05 exported.
## Differences from the original chunk, both deliberate:
##  * the species list is taken from THIS soil's metadata. The original read a
##    `species_list` left over from an earlier chunk (its definition was commented out),
##    which is how the mixed-system bacterial export ended up without the control.
##  * the sign convention is written into the file: logFC = log2(neighbour / control),
##    and `comparison` keeps the original "wheat.x.wheat vs <neighbour>" string so the
##    existing readers still parse it.
suppressPackageStartupMessages({library(edgeR); library(phyloseq)})

edger_vs_control <- function(ps, soil_code, control = "wheat x wheat", min_samples = 5,
                             p_export = 0.05) {
  sub  <- prune_samples(sample_data(ps)$soil == soil_code, ps)
  sub  <- filter_taxa(sub, function(x) sum(x > 0) >= min_samples, TRUE)
  meta <- data.frame(sample_data(sub), check.names = FALSE)
  counts <- round(as(otu_table(sub), "matrix"))
  if (!taxa_are_rows(sub)) counts <- t(counts)
  grp  <- factor(make.names(as.character(meta$species)))
  ctrl <- make.names(control)
  stopifnot(ctrl %in% levels(grp))
  grp  <- relevel(grp, ref = ctrl)
  dge  <- calcNormFactors(DGEList(counts = counts, group = grp), method = "TMM")
  design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
  dge  <- estimateDisp(dge, design)
  fit  <- glmQLFit(dge, design)
  out  <- lapply(setdiff(levels(grp), ctrl), function(sp) {
    con <- makeContrasts(contrasts = paste0("`", sp, "` - `", ctrl, "`"), levels = design)
    tt  <- topTags(glmQLFTest(fit, contrast = con), n = Inf)$table
    tt  <- tt[tt$PValue < p_export, , drop = FALSE]
    if (!nrow(tt)) return(NULL)
    tt$taxa <- rownames(tt); tt$comparison <- paste(ctrl, "vs", sp); tt
  })
  res <- do.call(rbind, out)
  fix <- c("Capsella.bursa.pastoris" = "Capsella.bursa-pastoris",
           "Legousia.speculum.veneris" = "Legousia.speculum-veneris",
           "Trigonella.foenum.graecum" = "Trigonella.foenum-graecum")
  for (k in names(fix)) res$comparison <- gsub(k, fix[[k]], res$comparison, fixed = TRUE)
  attr(res, "n_taxa_tested") <- nrow(counts)
  res
}
