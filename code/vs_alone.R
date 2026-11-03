## Wheat-alone check for Fig. 2 (Methods S3 (iii), Table S3 sheet vs_wheat_alone), called from Paper_figures.ipynb 3.1.2 via control_size.py.
## Each weed against wheat alone with the Fig. 2 rule: >= 2 of DESeq2 (BH), edgeR (BH), metagenomeSeq (Bonferroni, as in the paper's
## exports); adjusted P < 0.05; log2FC > 2. Wheat alone is the largest focal plant: a Fig. 2 ASV that also passes here is enriched relative
## to both a smaller (wheat x wheat) and a larger focal plant, and its enrichment cannot reflect plant size.
## Usage: Rscript vs_alone.R <matrix_dir> <out_dir>   (c. 10 min)
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(metagenomeSeq)})
args <- commandArgs(TRUE); MATRIX_DIR <- args[1]; OUT <- args[2]; dir.create(OUT, showWarnings = FALSE)
SYN <- c("Centaurea cyanus" = "Cyanus segetum", "Festuca arundinacea" = "Lolium arundinaceum", "Bromus sterilis" = "Anisantha sterilis")
build <- function(k, n) {
  otu  <- read.delim(file.path(MATRIX_DIR, sprintf("otu_%s_table_R_%dASVs.csv", k, n)), row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(MATRIX_DIR, sprintf("metadata_%s_%dASVs.csv", k, n)), check.names = FALSE)
  rownames(meta) <- meta$sample; keep <- intersect(colnames(otu), rownames(meta)); meta <- meta[keep, ]
  hit <- meta$species %in% names(SYN); meta$species[hit] <- SYN[meta$species[hit]]
  list(counts = round(as.matrix(otu[, keep])), meta = meta)
}
SOIL <- c(resistant = "Cereal-dominated (AF193)", Permissive = "Mixed-system (AF190)")
rows <- list()
for (kk in list(c("bacteria", 1314), c("fungi", 647))) {
  d <- build(kk[1], as.integer(kk[2]))
  for (sc in names(SOIL)) {
    m <- d$meta[d$meta$soil == sc, ]; X <- d$counts[, rownames(m)]; X <- X[rowSums(X > 0) >= 5, ]
    grp <- relevel(factor(make.names(m$species)), ref = "wheat.alone"); lv <- levels(grp)
    weeds <- setdiff(lv, c("wheat.alone", "wheat.x.wheat"))
    dds <- DESeqDataSetFromMatrix(X, data.frame(g = grp), ~ g)
    dds <- estimateSizeFactors(dds, type = "poscounts"); dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = TRUE)
    dge <- edgeR::calcNormFactors(DGEList(X, group = grp), method = "TMM"); des <- model.matrix(~0 + grp); colnames(des) <- lv
    dge <- estimateDisp(dge, des); fit <- glmQLFit(dge, des)
    for (w in weeds) {
      r1 <- as.data.frame(results(dds, contrast = c("g", w, "wheat.alone")))
      cv <- rep(0, length(lv)); cv[lv == w] <- 1; cv[lv == "wheat.alone"] <- -1
      r2 <- topTags(glmQLFTest(fit, contrast = cv), n = Inf, sort.by = "none")$table
      s  <- m$species %in% c(m$species[make.names(m$species) == w][1], "wheat alone"); Xs <- X[, rownames(m)[s]]; Xs <- Xs[rowSums(Xs > 0) >= 2, ]
      gs <- relevel(factor(ifelse(make.names(m$species[s]) == w, "weed", "alone")), ref = "alone")
      mr <- newMRexperiment(Xs); mr <- cumNorm(mr, p = cumNormStatFast(mr))
      fz <- try(fitZig(mr, model.matrix(~ gs), useCSSoffset = TRUE, control = zigControl(maxit = 30, verbose = FALSE)), silent = TRUE)
      r3 <- if (inherits(fz, "try-error")) NULL else MRcoefs(fz, coef = "gsweed", adjustMethod = "bonferroni", number = nrow(Xs))
      a <- rownames(X)
      L3 <- if (is.null(r3)) rep(NA, length(a)) else r3[[1]][match(a, rownames(r3))]
      Q3 <- if (is.null(r3)) rep(NA, length(a)) else r3$adjPvalues[match(a, rownames(r3))]
      df <- data.frame(kingdom = kk[1], soil = SOIL[[sc]], neighbour = m$species[make.names(m$species) == w][1], ASV = a,
                       DESeq2_l = r1$log2FoldChange, DESeq2_q = r1$padj, edgeR_l = r2$logFC, edgeR_q = r2$FDR, mgs_l = L3, mgs_q = Q3)
      ok <- function(l, q) !is.na(l) & !is.na(q) & l > 2 & q < 0.05
      df$n_methods <- ok(df$DESeq2_l, df$DESeq2_q) + ok(df$edgeR_l, df$edgeR_q) + ok(df$mgs_l, df$mgs_q)
      rows[[length(rows) + 1]] <- df
    }
    cat(kk[1], sc, "done\n")
  }
}
R <- do.call(rbind, rows); write.csv(R, file.path(OUT, "vs_wheat_alone_all.csv"), row.names = FALSE)
cat("calls (>=2 methods):", sum(R$n_methods >= 2), "\n")
