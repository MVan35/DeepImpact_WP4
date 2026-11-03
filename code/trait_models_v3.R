## =============================================================================
## trait_models_v3.R  --  Table 2, Notes S3 (traits) and the Mantel test
## -----------------------------------------------------------------------------
## Usage (from this folder):
##   Rscript trait_models_v3.R <data_dir> <trait_file> <tree_file> <out_dir> [nperm=999] [cores=4]
##   tree_file : weeds + wheat tree (V.PhyloMaker2, GBOTB.extended.WP, scenario S3), supplied as
##               weeds_wheat_VPM2_WP_S3.tre; rebuilt with V.PhyloMaker2 when that package is installed.
##
## Per kingdom x soil (weed pots only; the controls carry no weed traits):
##   1. Composition: one db-RDA per species-level trait on Bray-Curtis (Lingoes correction, as
##      vegan::dbrda(add = TRUE)); pseudo-F; P from reassigning trait values among the 28 species
##      (each species keeps its pots). Distance to wheat takes four values and is also tested as a
##      two-level contrast (monocots: Poaceae + Juncus vs eudicots).
##   2. Joint model of the seven traits: constrained R2, and its P from reassigning whole trait
##      profiles among species; species R2 (species as a factor) for comparison.
##   3. Richness: species means of log richness on the standardised traits (one model per cell).
##   4. Mantel test: weed x weed cophenetic distance vs species-centroid Bray-Curtis (9,999 perm.).
##   5. VIF of the traits at species level.
## =============================================================================
suppressPackageStartupMessages({library(vegan); library(parallel); library(ape)})
args   <- commandArgs(trailingOnly = TRUE)
DATA   <- args[1]; TRAITF <- args[2]; TREEF <- args[3]; OUT <- args[4]
NPERM  <- as.integer(if (length(args) >= 5) args[5] else 999)
CORES  <- as.integer(if (length(args) >= 6) args[6] else 4)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20261016)

SOILS <- c(resistant = "Cereal-dominated (AF193)", Permissive = "Mixed-system (AF190)")
JOINT  <- c("Distance", "Segetal", "Nitrogen", "Light", "SRL", "RDMC", "RootShoot")
TRAITS <- c(JOINT, "Monocot")

## ---- 1. phylogeny ------------------------------------------------------------
tr <- read.csv(TRAITF, check.names = FALSE)
meta0 <- read.csv(file.path(DATA, "metadata_bacteria_1314ASVs.csv"))
fam <- unique(meta0[!meta0$species %in% c("wheat alone", "wheat x wheat"), c("species", "family")])
splist <- rbind(data.frame(species = fam$species, genus = sub(" .*", "", fam$species), family = fam$family),
                data.frame(species = "Triticum aestivum", genus = "Triticum", family = "Poaceae"))
if (requireNamespace("V.PhyloMaker2", quietly = TRUE)) {
  phy <- V.PhyloMaker2::phylo.maker(splist, tree = V.PhyloMaker2::GBOTB.extended.WP,
                                    nodes = V.PhyloMaker2::nodes.info.1.WP, scenarios = "S3")$scenario.3
  write.tree(phy, file.path(OUT, "weeds_wheat_VPM2_WP_S3.tre"))
} else phy <- read.tree(TREEF)
coph <- cophenetic(phy); rownames(coph) <- colnames(coph) <- gsub("_", " ", rownames(coph))
tr$Distance <- coph[tr$species, "Triticum aestivum"]
tr$Monocot  <- as.numeric(tr$Distance < 250)          # Poaceae and Juncus (<= 188 Myr) vs eudicots (271.5 Myr)
dist_check <- data.frame(species = tr$species, published = tr$phy_dist_old, recomputed = tr$Distance,
                         changed = abs(tr$phy_dist_old - tr$Distance) > 1e-3)
write.csv(dist_check, file.path(OUT, "check_phylogenetic_distance.csv"), row.names = FALSE)
message("Distance changed for: ", paste(dist_check$species[dist_check$changed], collapse = ", "))

## ---- helpers: db-RDA by linear algebra (identical F to vegan::dbrda) ----------
gower_lingoes <- function(D) {
  D <- as.matrix(D); n <- nrow(D); J <- diag(n) - 1 / n
  G <- J %*% (-0.5 * D^2) %*% J
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -1e-10 * max(abs(ev))) {               # Lingoes (dbrda add = TRUE)
    D2 <- sqrt(D^2 + 2 * (-min(ev))); diag(D2) <- 0
    G <- J %*% (-0.5 * D2^2) %*% J
  }
  G
}
ssq <- function(X, G) { Q <- qr.Q(qr(cbind(1, X))); sum(Q * (G %*% Q)) }   # trace(H G)
pseudoF <- function(x, G) {
  n <- nrow(G); s <- ssq(x, G); tot <- sum(diag(G)); p <- NCOL(x)
  c(F = (s / p) / ((tot - s) / (n - p - 1)), R2 = s / tot)
}

## ---- per cell ------------------------------------------------------------------
run_cell <- function(kg, soil) {
  tag  <- if (kg == "bacteria") "1314" else "647"
  otu  <- read.delim(file.path(DATA, sprintf("otu_%s_table_R_%sASVs.csv", kg, tag)), row.names = 1, check.names = FALSE)
  meta <- read.csv(file.path(DATA, sprintf("metadata_%s_%sASVs.csv", kg, tag)))
  rownames(meta) <- meta$sample
  X <- t(as.matrix(otu)); meta <- meta[intersect(rownames(X), rownames(meta)), ]
  m <- meta[meta$soil == soil & !meta$species %in% c("wheat alone", "wheat x wheat"), ]
  X <- X[rownames(m), , drop = FALSE]
  T <- tr[match(m$species, tr$species), ]
  dat <- data.frame(species = m$species, Btot = m$Adventice_dry_mass,
                    Distance = T$Distance, Segetal = T$Segetal, Nitrogen = T$EllenN, Light = T$EllenL,
                    SRL = T$SRL, RDMC = T$RDMC, RootShoot = T$RootShoot, Monocot = T$Monocot,
                    rich = rowSums(X > 0), row.names = rownames(m))
  ## composition, one trait at a time
  comp <- do.call(rbind, lapply(TRAITS, function(v) {
    ok <- !is.na(dat[[v]]); d <- dat[ok, ]; G <- gower_lingoes(vegdist(X[ok, ], "bray"))
    obs <- pseudoF(d[[v]], G)
    ## pot-level permutation (free), as vegan::anova.cca does for a single term
    Fpot <- vapply(seq_len(NPERM), function(i) pseudoF(sample(d[[v]]), G)["F"], 0)
    ## species-level permutation: trait values reassigned among species
    spv <- tapply(d[[v]], d$species, `[`, 1); sps <- names(spv)
    Fsp <- vapply(seq_len(NPERM), function(i) {
      perm <- setNames(sample(spv), sps); pseudoF(perm[d$species], G)["F"] }, 0)
    Psp <- (sum(Fsp >= obs["F"]) + 1) / (NPERM + 1); nsp <- length(sps)
    data.frame(kingdom = kg, soil = SOILS[soil], trait = v, n_pots = nrow(d), n_species = nsp,
               F = obs["F"], R2 = obs["R2"], P_pot = (sum(Fpot >= obs["F"]) + 1) / (NPERM + 1), P_species = Psp)
  }))
  ## joint model of the seven traits; null = whole trait profiles reassigned among species
  cc <- complete.cases(dat[, JOINT]); dj <- dat[cc, ]; Gj <- gower_lingoes(vegdist(X[cc, ], "bray"))
  Z  <- scale(dj[, JOINT]); joint <- pseudoF(Z, Gj)
  prof <- unique(cbind(species = dj$species, as.data.frame(Z))); rownames(prof) <- prof$species
  R2null <- vapply(seq_len(NPERM), function(i) {
    s <- setNames(sample(prof$species), prof$species)
    pseudoF(as.matrix(prof[s[dj$species], JOINT]), Gj)["R2"] }, 0)
  R2sp <- pseudoF(model.matrix(~ factor(dj$species))[, -1], Gj)["R2"]
  jt <- data.frame(kingdom = kg, soil = SOILS[soil], n_pots = sum(cc), n_species = length(unique(dj$species)),
                   R2_total = joint["R2"], F_total = joint["F"], R2_null_mean = mean(R2null),
                   P_species = (sum(R2null >= joint["R2"]) + 1) / (NPERM + 1), R2_species = R2sp,
                   expected_share = length(JOINT) / (length(unique(dj$species)) - 1))
  ## vegan cross-check of one F value (guards the linear algebra above)
  ok <- !is.na(dat$Distance)
  vf <- anova(dbrda(vegdist(X[ok, ], "bray") ~ Distance, data = dat[ok, ], add = TRUE), permutations = 0)$F[1]
  stopifnot(abs(vf - comp$F[comp$trait == "Distance"]) < 1e-6)
  ## richness on species means (the pot-level model would treat pots as replicates of species traits)
  spm <- aggregate(dat[, c("rich", JOINT)], by = list(species = dat$species), FUN = mean)
  spm <- spm[complete.cases(spm), ]; spm$lr <- log(spm$rich); spm[, JOINT] <- scale(spm[, JOINT])
  sfit <- summary(lm(as.formula(paste("lr ~", paste(JOINT, collapse = " + "))), data = spm))$coefficients
  rich <- data.frame(kingdom = kg, soil = SOILS[soil], trait = JOINT, n_species = nrow(spm),
                     coef_std = sfit[JOINT, 1], P = sfit[JOINT, 4])
  ## Mantel: weed x weed phylogeny vs centroid dissimilarity
  cen <- rowsum(X / rowSums(X), m$species); cen <- cen / as.vector(table(m$species)[rownames(cen)])
  sp <- rownames(cen); Dm <- vegdist(cen, "bray"); Dp <- as.dist(coph[sp, sp])
  mt <- mantel(Dm, Dp, permutations = 9999)
  mant <- data.frame(kingdom = kg, soil = SOILS[soil], n_species = length(sp), mantel_r = mt$statistic, P = mt$signif)
  list(comp = comp, joint = jt, rich = rich, mantel = mant)
}
cells <- expand.grid(kg = c("bacteria", "fungi"), soil = names(SOILS), stringsAsFactors = FALSE)
R <- mclapply(seq_len(nrow(cells)), function(i) { set.seed(20261016 + i)   # one seed per cell: reproducible with any number of cores
  run_cell(cells$kg[i], cells$soil[i]) }, mc.cores = CORES)
for (k in c("comp", "joint", "rich", "mantel"))
  write.csv(do.call(rbind, lapply(R, `[[`, k)), file.path(OUT, sprintf("trait_%s.csv", k)), row.names = FALSE)
## collinearity at species level
sp <- data.frame(Distance = tr$Distance, Segetal = tr$Segetal, Nitrogen = tr$EllenN, Light = tr$EllenL,
                 SRL = tr$SRL, RDMC = tr$RDMC, RootShoot = tr$RootShoot); sp <- sp[complete.cases(sp), ]
vif <- sapply(JOINT, function(v) 1 / (1 - summary(lm(as.formula(paste(v, "~ .")), data = sp[, JOINT]))$r.squared))
write.csv(data.frame(trait = names(vif), VIF = vif), file.path(OUT, "trait_vif.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
message("done")
