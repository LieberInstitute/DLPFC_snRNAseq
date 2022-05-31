library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")

## Choose Correction 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
correction = args[1]

## Load empty-free sce data
load(here("processed-data", "03_build_sce","sce_uncorrected.Rdata"), verbose = TRUE)
sce <- sce_uncorrected

dim(sce)
# [1] 36601 77604

stopifnot(correction %in% colnames(colData(sce)))

#### Preform Batch Correction with MNN ####
table(sce[[correction]])

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(519)

## what about merge order? maybe use auto.merge=TRUE
message("running fast MNN - ", Sys.time())

mnn <- fastMNN(sce, batch=sce[[correction]],
                     auto.merge=TRUE,
                     subset.row=rowData(sce)$hvgs, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())

message("Done MNN  ", Sys.time())

# Add them to the SCE, as well as the metadata
reducedDim(sce, "PCA") <- reducedDim(mnn, "corrected") # 100 components 
metadata(sce) <- metadata(mnn)

## Check out lost variance (looking for < 10%)
100 * metadata(mnn)$merge.info$lost.var

## Should we find optimal PC space? - use all 100 for now

set.seed(109)

message("running TSNE - ", Sys.time())
sce <- runTSNE(sce, dimred = pca)

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = pca)

message("Done TSNE + UMAP - Saving data...", Sys.time())
  
save(sce, file = here("processed-data", "03_build_sce",paste0("sce_MNN_",correction,".Rdata")))

# sgejobs::job_single('normalize_MNN_round', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript normalize_MNN.R round")
# sgejobs::job_single('normalize_MNN_subject', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript normalize_MNN.R subject")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
