## Order-independence of the Table 1 partition (Methods S2).
## Usage: Rscript permanova_order.R <data_dir> <out_csv>
suppressPackageStartupMessages(library(vegan))
a <- commandArgs(TRUE); set.seed(1)
res <- list()
for (kg in c("bacteria", "fungi")) {
  tag <- if (kg == "bacteria") "1314" else "647"
  X <- t(as.matrix(read.delim(file.path(a[1], sprintf("otu_%s_table_R_%sASVs.csv", kg, tag)), row.names = 1, check.names = FALSE)))
  m <- read.csv(file.path(a[1], sprintf("metadata_%s_%sASVs.csv", kg, tag))); rownames(m) <- m$sample
  m <- m[rownames(X), ]; D <- vegdist(X, "bray")
  f <- list(published_order = D ~ soil * species, reversed_order = D ~ species * soil)
  for (nm in names(f)) { r <- adonis2(f[[nm]], data = m, permutations = 999, by = "terms")
    res[[length(res) + 1]] <- data.frame(kingdom = kg, model = nm, term = rownames(r), df = r$Df, R2 = r$R2, F = r$F, P = r[["Pr(>F)"]]) }
  r <- adonis2(D ~ soil + species, data = m, permutations = 999, by = "margin")
  res[[length(res) + 1]] <- data.frame(kingdom = kg, model = "marginal_main_effects", term = rownames(r), df = r$Df, R2 = r$R2, F = r$F, P = r[["Pr(>F)"]])
}
write.csv(do.call(rbind, res), a[2], row.names = FALSE)
