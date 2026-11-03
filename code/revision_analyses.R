## =============================================================================
## revision_analyses.R  --  analyses added after the two statistical audits
## -----------------------------------------------------------------------------
## Usage (from this folder):
##   Rscript revision_analyses.R <data_dir> <out_dir> [nperm=9999]
##
##   Table 1   within-soil neighbour effect: adonis2(Bray-Curtis ~ species) per kingdom x soil
##   Fig. 4    direct microbiota-growth test: db-RDA of species-centroid composition on species-mean
##             wheat biomass, 28 weeds per kingdom x soil, with conditional versions
##             (monocot/eudicot; monocot + legume; neighbour biomass) and without legumes;
##             site scores for the figure
##             and a per-taxon Spearman screen at the same unit (BH within kingdom x soil)
##   Fig. S6   Glomeromycota share by soil (Wilcoxon, pots) and Fabaceae vs other weeds
##             (Wilcoxon on species means, per soil and pooled)
## =============================================================================
suppressPackageStartupMessages({library(vegan)})
args <- commandArgs(trailingOnly = TRUE)
DATA <- args[1]; OUT <- args[2]; NPERM <- as.integer(if (length(args) >= 3) args[3] else 9999)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); set.seed(20261102)
SOILS <- c(resistant = "Cereal-dominated (AF193)", Permissive = "Mixed-system (AF190)")
CTRL  <- c("wheat alone", "wheat x wheat")
SYN   <- c("Centaurea cyanus" = "Cyanus segetum", "Festuca arundinacea" = "Lolium arundinaceum", "Bromus sterilis" = "Anisantha sterilis")
LEGUMES <- c("Lotus corniculatus", "Medicago lupulina", "Trigonella foenum-graecum", "Vicia cracca")
load_kingdom <- function(kg) {
  tag  <- if (kg == "bacteria") "1314" else "647"
  otu  <- read.delim(file.path(DATA, sprintf("otu_%s_table_R_%sASVs.csv", kg, tag)), row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(DATA, sprintf("metadata_%s_%sASVs.csv", kg, tag)), check.names = FALSE)
  rownames(meta) <- meta$sample; hit <- meta$species %in% names(SYN); meta$species[hit] <- SYN[meta$species[hit]]
  X <- t(as.matrix(otu)); keep <- intersect(rownames(X), rownames(meta))
  list(X = X[keep, ], meta = meta[keep, ])
}
K <- list(bacteria = load_kingdom("bacteria"), fungi = load_kingdom("fungi"))

## ---- Table 1: within-soil neighbour effect -------------------------------------
within <- do.call(rbind, lapply(names(K), function(kg) do.call(rbind, lapply(names(SOILS), function(sc) {
  m <- K[[kg]]$meta[K[[kg]]$meta$soil == sc, ]; X <- K[[kg]]$X[rownames(m), ]
  do.call(rbind, lapply(c("all treatments", "weeds only"), function(set) {
    mm <- if (set == "weeds only") m[!m$species %in% CTRL, ] else m
    a <- adonis2(vegdist(X[rownames(mm), ], "bray") ~ species, data = mm, permutations = 999)
    data.frame(kingdom = kg, soil = SOILS[[sc]], set = set, n_pots = nrow(mm), df = a$Df[1],
               R2_pct = 100 * a$R2[1], F = a$F[1], P = a$`Pr(>F)`[1])
  }))
}))))
write.csv(within, file.path(OUT, "table1_within_soil.csv"), row.names = FALSE)

## ---- Fig. 4: direct microbiota-growth test ---------------------------------------
direct <- list(); scores <- list(); taxa <- list()
for (kg in names(K)) for (sc in names(SOILS)) {
  m <- K[[kg]]$meta[K[[kg]]$meta$soil == sc & !is.na(K[[kg]]$meta$Dry_wheat_mass) & !K[[kg]]$meta$species %in% CTRL, ]
  X <- K[[kg]]$X[rownames(m), ]; X <- X / rowSums(X)
  cen <- rowsum(X, m$species) / as.vector(table(m$species)[rownames(rowsum(X, m$species))])
  sp  <- rownames(cen)
  d <- data.frame(species = sp,
                  biomass = tapply(m$Dry_wheat_mass, m$species, mean)[sp],
                  Btot    = tapply(m$Adventice_dry_mass, m$species, mean, na.rm = TRUE)[sp],
                  family  = tapply(m$family, m$species, `[`, 1)[sp])
  d$Monocot <- as.numeric(d$family %in% c("Poaceae", "Juncaceae")); d$Legume <- as.numeric(d$species %in% LEGUMES)
  D <- vegdist(cen, "bray")
  fit <- function(f, dd = d, DD = D) { a <- anova(dbrda(f, data = dd), permutations = how(nperm = NPERM)); c(F = a$F[1], P = a$`Pr(>F)`[1]) }
  m0 <- dbrda(D ~ biomass, data = d); r2 <- RsquareAdj(m0)$r.squared
  u  <- fit(D ~ biomass)
  c1 <- fit(D ~ biomass + Condition(Monocot))
  c2 <- fit(D ~ biomass + Condition(Monocot + Legume))
  okb <- !is.na(d$Btot); Db <- vegdist(cen[okb, ], "bray")
  c3 <- fit(Db ~ biomass + Condition(Btot), d[okb, ], Db)
  nl <- d$Legume == 0; Dn <- vegdist(cen[nl, ], "bray"); c4 <- fit(Dn ~ biomass, d[nl, ], Dn)
  direct[[length(direct) + 1]] <- data.frame(kingdom = kg, soil = SOILS[[sc]], n_species = nrow(d), R2_pct = 100 * r2,
    F = u["F"], P = u["P"], P_given_monocot = c1["P"], P_given_monocot_legume = c2["P"], P_given_Btot = c3["P"],
    R2_adj_pct = 100 * RsquareAdj(m0)$adj.r.squared,
    R2_without_legumes_pct = 100 * RsquareAdj(dbrda(Dn ~ biomass, data = d[nl, ]))$r.squared,
    P_without_legumes = c4["P"], r_biomass_monocot = cor(d$biomass, d$Monocot))
  ## per-taxon screen at the same unit: Spearman of species-mean relative abundance with species-mean biomass
  elig <- colnames(X)[colMeans(X > 0) >= 0.10 & colMeans(X) >= 1e-4]
  st <- t(sapply(elig, function(a) { ct <- suppressWarnings(cor.test(cen[, a], d$biomass, method = "spearman", exact = FALSE)); c(ct$estimate, ct$p.value) }))
  taxa[[length(taxa) + 1]] <- data.frame(kingdom = kg, soil = SOILS[[sc]], ASV = elig, rho = st[, 1], P = st[, 2], q_BH = p.adjust(st[, 2], "BH"))
  sc1 <- scores(m0, display = "sites", choices = 1:2, scaling = 1)
  scores[[length(scores) + 1]] <- data.frame(kingdom = kg, soil = SOILS[[sc]], species = sp, family = d$family,
    biomass = d$biomass, CAP1 = sc1[, 1], MDS1 = sc1[, 2],
    eig_CAP1_pct = 100 * m0$CCA$eig[1] / m0$tot.chi, eig_MDS1_pct = 100 * m0$CA$eig[1] / m0$tot.chi)
}
direct <- do.call(rbind, direct); direct$q_BH <- p.adjust(direct$P, "BH")
taxa <- do.call(rbind, taxa); write.csv(taxa, file.path(OUT, "fig4_taxon_screen.csv"), row.names = FALSE)
write.csv(direct, file.path(OUT, "fig4_direct_test.csv"), row.names = FALSE)
write.csv(do.call(rbind, scores), file.path(OUT, "fig4_scores.csv"), row.names = FALSE)

## ---- Fig. S6: Glomeromycota --------------------------------------------------------
tax <- read.delim(file.path(DATA, "tax_fungi_R_647ASVs.csv"), row.names = 1, check.names = FALSE)
glom <- rownames(tax)[grepl("^Glomeromycota", tax$Phylum)]
m <- K$fungi$meta; X <- K$fungi$X[rownames(m), ]; m$glom <- 100 * rowSums(X[, intersect(glom, colnames(X))]) / rowSums(X)
w <- wilcox.test(glom ~ soil, data = m)
spm <- aggregate(glom ~ soil + species + family, data = m[!m$species %in% CTRL, ], FUN = mean)
fab <- function(s) { x <- if (is.null(s)) spm else spm[spm$soil == s, ]; wilcox.test(x$glom[x$family == "Fabaceae"], x$glom[x$family != "Fabaceae"])$p.value }
glo <- data.frame(test = c("soil (pots)", "Fabaceae vs other weeds, species means, cereal-dominated", "same, mixed-system", "same, both soils pooled"),
  stat = c(sprintf("mean %.2f vs %.2f %%; median %.2f vs %.2f %% (mixed vs cereal)",
                   mean(m$glom[m$soil == "Permissive"]), mean(m$glom[m$soil == "resistant"]),
                   median(m$glom[m$soil == "Permissive"]), median(m$glom[m$soil == "resistant"])), "", "", ""),
  n = c(nrow(m), sum(spm$soil == "resistant"), sum(spm$soil == "Permissive"), nrow(spm)),
  P = c(w$p.value, fab("resistant"), fab("Permissive"), fab(NULL)))
write.csv(glo, file.path(OUT, "figS6_glomeromycota_tests.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_revision.txt"))
message("done")
