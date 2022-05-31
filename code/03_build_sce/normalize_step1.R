library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")

## Load empty-free sce data
load(here("processed-data", "sce", "sce_no_empty_droplets.Rdata"), verbose = TRUE)

## Drop auto-drop nuclei
table(sce$discard_auto)
# FALSE  TRUE 
# 77604  7152

sce <- sce[,!sce$discard_auto]
dim(sce)
# [1] 36601 77604

#### Rescale ####
# Use `multiBatchNorm()` to compute log-normalized counts
sce <- multiBatchNorm(sce, batch=sce$round)

# Find HVGs
geneVar <- modelGeneVar(sce)
chosen.hvgs <- geneVar$bio > 0
sum(chosen.hvgs)
# [1] 13254

rowData(sce)$hvgs <- chosen.hvgs

#### How does the un-normalized data look? ####
sce_uncorrected <- runPCA(sce, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())

sce_uncorrected <- runTSNE(uncorrected, dimred="PCA")

save(sce_uncorrected, file = here("processed-data", "03_build_sce","sce_uncorrected.Rdata"))

# sgejobs::job_single('normalize_step1', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize_step1.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
