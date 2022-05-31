library("SingleCellExperiment")
library("harmony")
library("here")
library("sessioninfo")

## Choose Correction 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
correction = args[1]

stopifnot(correction %in% colnames(colData(sce)))

## Load empty-free sce data
load(here("processed-data", "03_build_sce","sce_uncorrected.Rdata"), verbose = TRUE)

message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce_uncorrected, correction, lambda = .1, verbose = TRUE)

set.seed(109)

message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = 'HARMONY')

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = 'HARMONY')

message("Done TSNE + UMAP - Saving data...", Sys.time())

save(sce, file = here("processed-data", "03_build_sce",paste0("sce_harmony_",correction,".Rdata")))


# sgejobs::job_single('normalize_harmony_round', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript normalize_harmony.R round")
# sgejobs::job_single('normalize_harmony_subject', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript normalize_harmony.R subject")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
