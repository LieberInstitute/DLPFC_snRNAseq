library("SingleCellExperiment")
library("zellkonverter")
library("here")
library("sessioninfo")

load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

annData <- SCE2AnnData(sce)

