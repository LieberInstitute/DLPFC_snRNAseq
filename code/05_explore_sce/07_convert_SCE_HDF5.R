library("SingleCellExperiment")
library("here")
library("sessioninfo")

## load data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

counts(sce)
# names(3): counts binomial_deviance_residuals logcounts

# Just save log counts
assays(sce)$counts <- NULL
assays(sce)$binomial_deviance_residuals <- NULL

message("Object Size: ", lobstr::obj_size(sce))
# 181.85 MB

## print details 
message("\n **** SCE details ****")
sce

message("\n **** colData ****")
colData(sce)

message("\n **** rowData ****")
rowData(sce)

saveHDF5SummarizedExperiment(sce, here("processed-data", "sce", "sce_DLPFC_annotated"), replace = TRUE)

# sgejobs::job_single('07_convert_SCE_HDF5', create_shell = TRUE, memory = '25G', command = "Rscript 07_convert_SCE_HDF5.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
