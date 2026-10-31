#!/usr/bin/env Rscript
# Reconcile pcrelate_to_mtx.py against GENESIS::pcrelate on the selftest fixture (2026-10-31).
#
# Two implementations of one estimator must agree pair for pair. Run the python selftest with
# --keep DIR (it writes the simulated .bed/.bim/.fam, the PCs, the unrelated set and EVERY pair's
# kinship), then this script on DIR/sim. GENESIS does not build against the pixi R (gdsfmt fails to
# compile, the Fortran wrapper is broken); install it as conda binaries instead:
#   conda create -p <env> -c conda-forge -c bioconda bioconductor-genesis bioconductor-snprelate \
#       bioconductor-gwastools bioconductor-gdsfmt r-base
#   <env>/bin/Rscript harness/tools/pcrelate_genesis_recon.R DIR/sim
# Measured 2026-10-31 on the 471 x 16,000 fixture: 110,685 pairs, max |diff| 7.2e-4, median
# 4.7e-5, correlation 0.99997; inbreeding max |diff| 0.02 (GENESIS's kinSelf is a different
# form, which is why the python tool's diagonal defaults to 1). GENESIS is given the same PCs,
# the same training set and maf.thresh 0.01, small.samp.correct FALSE.
suppressPackageStartupMessages({library(SNPRelate); library(GENESIS); library(GWASTools); library(gdsfmt)})
args <- commandArgs(TRUE)
if (length(args) != 1) stop("usage: pcrelate_genesis_recon.R DIR/sim  (the --keep prefix of the python selftest)")
pre <- args[1]
gds <- paste0(pre, ".gds")
if (file.exists(gds)) unlink(gds)
snpgdsBED2GDS(paste0(pre, ".bed"), paste0(pre, ".fam"), paste0(pre, ".bim"), gds, cvt.chr = "int", verbose = FALSE)
pcs <- read.table(paste0(pre, ".pcs.tsv"), header = TRUE, comment.char = "", check.names = FALSE)
names(pcs)[1:2] <- c("FID", "IID")
pcm <- as.matrix(pcs[, grep("^PC", names(pcs))])
rownames(pcm) <- pcs$IID
unrel <- readLines(paste0(pre, ".unrel.txt"))
gd <- GenotypeData(GdsGenotypeReader(gds))
pcr <- pcrelate(GenotypeBlockIterator(gd, snpBlock = 4000), pcs = pcm, training.set = unrel,
                maf.thresh = 0.01, small.samp.correct = FALSE, verbose = FALSE)
kb <- pcr$kinBtwn
ib <- pcr$kinSelf
mine <- read.table(paste0(pre, ".pairs_all.tsv"), header = TRUE, colClasses = c("character", "character", "numeric"))
key <- function(a, b) ifelse(a > b, paste(a, b), paste(b, a))
m <- merge(data.frame(k = key(kb$ID1, kb$ID2), genesis = kb$kin), data.frame(k = key(mine$i, mine$j), mine = mine$kinship))
d <- m$mine - m$genesis
cat(sprintf("pairs compared %d (GENESIS %d, mine %d); max |diff| %.3g; median |diff| %.3g; cor %.6f\n",
            nrow(m), nrow(kb), nrow(mine), max(abs(d)), median(abs(d)), cor(m$mine, m$genesis)))
lab <- read.table(paste0(pre, ".labels.tsv"), header = TRUE, colClasses = c("character", "character", "numeric"))
f <- merge(data.frame(iid = ib$ID, genesis_f = ib$f), lab)
cat(sprintf("inbreeding: %d samples, max |diff| %.3g, inbred child GENESIS %.4f mine %.4f\n",
            nrow(f), max(abs(f$genesis_f - f$inbreeding)), f$genesis_f[f$label == "inbred"], f$inbreeding[f$label == "inbred"]))
ok <- nrow(m) == nrow(mine) && max(abs(d)) < 1e-3
cat(if (ok) "reconciliation OK\n" else "reconciliation FAILED\n")
quit(status = if (ok) 0 else 1)
