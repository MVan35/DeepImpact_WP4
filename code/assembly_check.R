## Turnover (betapart) and bacterial βNTI (iCAMP::bNTIn.p, from the NST authors) for Methods S2 / Fig. S7.
## Turnover runs anywhere betapart is installed; βNTI needs iCAMP (install.packages("iCAMP")) and the rpoB ASV tree.
## Usage: Rscript assembly_check.R <matrix_dir> <rpoB_tree.nwk> <out_dir> [rand = 999] [nworker = 4]
suppressPackageStartupMessages({library(betapart); library(ape)})
a <- commandArgs(TRUE); DF <- a[1]; TREE <- a[2]; OUT <- a[3]
RAND <- as.integer(if (length(a) >= 4) a[4] else 999); NW <- as.integer(if (length(a) >= 5) a[5] else 4)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(20261103)
SOIL <- c(resistant = "Cereal-dominated (AF193)", Permissive = "Mixed-system (AF190)")
load_k <- function(k, n) {
  X <- t(as.matrix(read.delim(file.path(DF, sprintf("otu_%s_table_R_%dASVs.csv", k, n)), row.names = 1, check.names = FALSE)))
  m <- read.csv(file.path(DF, sprintf("metadata_%s_%dASVs.csv", k, n))); rownames(m) <- m$sample
  m <- m[intersect(rownames(X), rownames(m)), ]; list(X = X[rownames(m), ], m = m)
}
## ---- 1. turnover share of Sorensen dissimilarity (betapart::beta.pair), per kingdom x soil
tv <- list()
for (kk in list(c("bacteria", 1314), c("fungi", 647))) {
  d <- load_k(kk[1], as.integer(kk[2]))
  for (s in names(SOIL)) {
    PA <- (d$X[d$m$soil == s, ] > 0) * 1; PA <- PA[, colSums(PA) > 0]
    bp <- beta.pair(PA, index.family = "sorensen")
    sim <- as.vector(bp$beta.sim); sor <- as.vector(bp$beta.sor)
    tv[[length(tv) + 1]] <- data.frame(kingdom = kk[1], soil = SOIL[[s]], pots = nrow(PA), mean_beta_sor = mean(sor), mean_beta_sim = mean(sim),
                                       mean_beta_sne = mean(as.vector(bp$beta.sne)), turnover_share_pct = 100 * mean(sim) / mean(sor),
                                       turnover_share_pairwise_pct = 100 * mean((sim / sor)[sor > 0]))
  }
}
tv <- do.call(rbind, tv); write.csv(tv, file.path(OUT, "turnover_betapart.csv"), row.names = FALSE); print(tv, digits = 4)
## ---- 2. bacterial βNTI per neighbour x soil group (iCAMP::bNTIn.p, abundance-weighted βMNTD, taxa shuffled across tips)
if (requireNamespace("iCAMP", quietly = TRUE) && file.exists(TREE)) {
  tr <- read.tree(TREE); d <- load_k("bacteria", 1314); keep <- intersect(colnames(d$X), tr$tip.label)
  tr <- keep.tip(tr, keep); pd <- cophenetic(tr); rows <- list()
  for (s in names(SOIL)) for (g in sort(unique(d$m$species[d$m$soil == s]))) {
    idx <- rownames(d$m)[d$m$soil == s & d$m$species == g]; if (length(idx) < 3) next
    cm <- d$X[idx, keep]; cm <- cm[, colSums(cm) > 0, drop = FALSE]
    b <- iCAMP::bNTIn.p(comm = cm, dis = pd[colnames(cm), colnames(cm)], nworker = NW, weighted = TRUE, rand = RAND)
    v <- as.matrix(b$index); v <- v[upper.tri(v)]
    rows[[length(rows) + 1]] <- data.frame(soil = SOIL[[s]], species = g, pots = length(idx), mean_bNTI = mean(v), sd_bNTI = sd(v),
                                           pct_pairs_above2 = 100 * mean(v > 2), pct_pairs_below_minus2 = 100 * mean(v < -2))
    cat(s, g, round(mean(v), 2), "\n")
  }
  bn <- do.call(rbind, rows); write.csv(bn, file.path(OUT, "bNTI_bacteria_iCAMP.csv"), row.names = FALSE)
  cat("groups with mean bNTI > 2:", sum(bn$mean_bNTI > 2), "of", nrow(bn), "\n")
} else message("iCAMP or tree not available: bNTI skipped")
