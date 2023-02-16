library("SingleCellExperiment")
library("here")
library("sessioninfo")

## load data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce", "sce_DLPFC"))
rhdf5::h5ls(here("processed-data", "sce", "sce_DLPFC","assays.h5"))
sce_h5 <- rhdf5::h5read(here("processed-data", "sce", "sce_DLPFC","assays.h5"), "assay001")
sce <- zellkonverter::readH5AD(here("processed-data", "sce", "sce_DLPFC","assays.h5"))

names(assays(sce))
# counts binomial_deviance_residuals logcounts

message("Size with all assasy: ")
lobstr::obj_size(sce)
# 311.00 MB

## print details
message("\n **** SCE details ****")
sce

message("\n **** colData ****")
colData(sce)

message("\n **** rowData ****")
rowData(sce)

message(Sys.time(), "- Saving Data")
saveHDF5SummarizedExperiment(sce, here("processed-data", "sce", "sce_DLPFC_annotated"), replace = TRUE)

# Just save log counts
assays(sce)$counts <- NULL
assays(sce)$binomial_deviance_residuals <- NULL

message("Size with only logcounts: ")
lobstr::obj_size(sce)
# 181.85 MB

message(Sys.time(), "- Saving Data")
saveHDF5SummarizedExperiment(sce, here("processed-data", "sce", "sce_DLPFC_annotated-logcounts_only"), replace = TRUE)
## to load
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce", "sce_DLPFC_annotated-logcounts_only"))
# sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce", "sce_DLPFC_annotated"))

# sgejobs::job_single('07_convert_SCE_HDF5', create_shell = TRUE, memory = '25G', command = "Rscript 07_convert_SCE_HDF5.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
