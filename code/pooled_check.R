## Control-size check for Fig. 2 (decision D3).
## (i)  Pooled contrast: all 28 weeds pooled vs wheat x wheat, per kingdom x soil, with the Fig. 2 rule
##      (>= 2 of DESeq2, edgeR, metagenomeSeq; BH q < 0.05; log2FC > 2). ASVs called here respond to
##      "any weed" (or to the larger focal wheat), not to a particular neighbour.
## (ii) Specificity contrast: for each Fig. 2 call (neighbour, ASV), that neighbour's pots vs the pots of the other
##      27 weeds pooled (DESeq2 and edgeR), one-sided, BH over the calls of that kingdom x soil.
##      This comparison does not involve the wheat x wheat control at all.
## Usage: Rscript pooled_check.R <matrix_dir> <fig2_values.csv> <out_dir>
suppressPackageStartupMessages({library(phyloseq); library(DESeq2); library(edgeR); library(metagenomeSeq)})
args <- commandArgs(TRUE); MATRIX_DIR <- args[1]; FIG2 <- args[2]; OUT <- args[3]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SYN <- c("Centaurea cyanus" = "Cyanus segetum", "Festuca arundinacea" = "Lolium arundinaceum", "Bromus sterilis" = "Anisantha sterilis")
CTRL <- c("wheat x wheat", "wheat alone")
build <- function(k, n) {
  otu  <- read.delim(file.path(MATRIX_DIR, sprintf("otu_%s_table_R_%dASVs.csv", k, n)), row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(MATRIX_DIR, sprintf("metadata_%s_%dASVs.csv", k, n)), check.names = FALSE)
  rownames(meta) <- meta$sample; keep <- intersect(colnames(otu), rownames(meta)); meta <- meta[keep, ]
  hit <- meta$species %in% names(SYN); meta$species[hit] <- SYN[meta$species[hit]]
  list(counts = round(as.matrix(otu[, keep])), meta = meta)
}
fig2 <- read.csv(FIG2, check.names = FALSE)
SOIL <- c(resistant = "Cereal-dominated (AF193)", Permissive = "Mixed-system (AF190)")
pooled_rows <- list(); spec_rows <- list()
for (kk in list(c("bacteria", 1314), c("fungi", 647))) {
  d <- build(kk[1], as.integer(kk[2]))
  for (sc in names(SOIL)) {
    idx <- d$meta$soil == sc; m <- d$meta[idx, ]; X <- d$counts[, rownames(m)]
    X <- X[rowSums(X > 0) >= 5, ]
    ## ---------- (i) pooled contrast
    grp <- factor(ifelse(m$species %in% CTRL, make.names(m$species), "weed"))
    grp <- relevel(grp, ref = "wheat.x.wheat")
    dds <- DESeqDataSetFromMatrix(X, data.frame(g = grp), ~ g)
    dds <- estimateSizeFactors(dds, type = "poscounts"); dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
    r1 <- results(dds, contrast = c("g", "weed", "wheat.x.wheat"))
    dge <- edgeR::calcNormFactors(DGEList(X, group = grp), method = "TMM"); des <- model.matrix(~0 + grp); colnames(des) <- levels(grp)
    dge <- estimateDisp(dge, des); fit <- glmQLFit(dge, des)
    r2 <- topTags(glmQLFTest(fit, contrast = makeContrasts(weed - wheat.x.wheat, levels = des)), n = Inf, sort.by = "none")$table
    mr <- newMRexperiment(X); mr <- cumNorm(mr, p = cumNormStatFast(mr))
    mod <- model.matrix(~ grp); fz <- fitZig(mr, mod, useCSSoffset = TRUE, control = zigControl(maxit = 50, verbose = FALSE))
    r3 <- MRcoefs(fz, coef = "grpweed", adjustMethod = "bonferroni", number = nrow(X)); cat("mgs cols:", colnames(r3), "\n")
    co <- function(x, fc, q) setNames(data.frame(x[[fc]], x[[q]], row.names = rownames(x)), c("l", "q"))
    T3 <- list(DESeq2 = co(as.data.frame(r1), "log2FoldChange", "padj"), edgeR = co(r2, "logFC", "FDR"),
               metagenomeSeq = setNames(data.frame(r3[[1]], r3$adjPvalues, row.names = rownames(r3)), c("l", "q")))
    call <- sapply(T3, function(t) { v <- rep(FALSE, nrow(X)); names(v) <- rownames(X)
                                     ok <- !is.na(t$q) & t$q < 0.05 & t$l > 2; v[rownames(t)[ok]] <- TRUE; v })
    pooled_rows[[paste(kk[1], sc)]] <- data.frame(kingdom = kk[1], soil = SOIL[[sc]], ASV = rownames(X),
        n_methods = rowSums(call), DESeq2_l = T3$DESeq2[rownames(X), "l"], DESeq2_q = T3$DESeq2[rownames(X), "q"],
        edgeR_l = T3$edgeR[rownames(X), "l"], edgeR_q = T3$edgeR[rownames(X), "q"],
        mgs_l = T3$metagenomeSeq[rownames(X), "l"], mgs_q = T3$metagenomeSeq[rownames(X), "q"])
    ## ---------- (ii) neighbour vs all other weeds pooled (two-level factor; no wheat x wheat)
    f2 <- fig2[fig2$kingdom == kk[1] & fig2$soil == SOIL[[sc]], ]
    if (!nrow(f2)) next
    w  <- !(m$species %in% CTRL); mw <- m[w, ]; Xw <- X[, rownames(mw)]; Xw <- Xw[rowSums(Xw > 0) >= 3, ]
    for (nb in unique(f2$neighbour)) {
      g2 <- factor(ifelse(mw$species == nb, "this", "other"), levels = c("other", "this"))
      dd <- DESeqDataSetFromMatrix(Xw, data.frame(g = g2), ~ g)
      dd <- estimateSizeFactors(dd, type = "poscounts"); dd <- DESeq(dd, test = "Wald", fitType = "local", quiet = TRUE)
      rd <- as.data.frame(results(dd, contrast = c("g", "this", "other")))
      dg <- edgeR::calcNormFactors(DGEList(Xw, group = g2), method = "TMM"); ds <- model.matrix(~ g2)
      dg <- estimateDisp(dg, ds); re <- glmQLFTest(glmQLFit(dg, ds), coef = 2)$table
      for (a in f2$ASV[f2$neighbour == nb]) {
        ok <- a %in% rownames(Xw)
        spec_rows[[length(spec_rows) + 1]] <- data.frame(kingdom = kk[1], soil = SOIL[[sc]], neighbour = nb, ASV = a,
          DESeq2_l = if (ok) rd[a, "log2FoldChange"] else NA, DESeq2_p = if (ok) rd[a, "pvalue"] else NA,
          edgeR_l = if (ok) re[a, "logFC"] else NA, edgeR_p = if (ok) re[a, "PValue"] else NA)
      }
    }
    cat(kk[1], sc, "done\n")
  }
}
P <- do.call(rbind, pooled_rows); S <- do.call(rbind, spec_rows)
for (k in unique(paste(S$kingdom, S$soil))) {
  i <- paste(S$kingdom, S$soil) == k
  S$DESeq2_q[i] <- p.adjust(ifelse(S$DESeq2_l[i] > 0, S$DESeq2_p[i] / 2, 1 - S$DESeq2_p[i] / 2), "BH")   # one-sided: higher under this neighbour
  S$edgeR_q[i]  <- p.adjust(ifelse(S$edgeR_l[i] > 0, S$edgeR_p[i] / 2, 1 - S$edgeR_p[i] / 2), "BH")
}
write.csv(P, file.path(OUT, "pooled_weeds_vs_wxw.csv"), row.names = FALSE)
write.csv(S, file.path(OUT, "fig2_neighbour_vs_other_weeds.csv"), row.names = FALSE)
cat("pooled ASVs called:", sum(P$n_methods >= 2), "| specificity rows:", nrow(S), "\n")
