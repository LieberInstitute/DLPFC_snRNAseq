library("SingleCellExperiment")
library("harmony")
library("scater")
library("here")
library("sessioninfo")

## Choose Correction 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
correction = args[1]

## Load empty-free sce data
load(here("processed-data", "03_build_sce","sce_uncorrected_glm.Rdata"), verbose = TRUE)
stopifnot(correction %in% colnames(colData(sce_uncorrected)))

## Run harmony
## needs PCA & GLMPCA? not sure which it acctauly uses? 
reducedDim(sce_uncorrected, "PCA") <- reducedDim(sce_uncorrected, "GLMPCA_approx")
# Error in RunHarmony.SingleCellExperiment(sce_uncorrected, reduction = "GLMPCA_approx",  : 
#                                            PCA must be computed before running Harmony.

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce_uncorrected, reduction="GLMPCA_approx", group.by.vars = correction, verbose = TRUE)

set.seed(602)

message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = 'HARMONY')

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = 'HARMONY')

message("Done TSNE + UMAP - Saving data...", Sys.time())

save(sce, file = here("processed-data", "03_build_sce",paste0("sce_harmony_",correction,".Rdata")))

## going with Sample
# sgejobs::job_single('05_harmony_correction_Sample', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 05_harmony_correction.R Sample")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
