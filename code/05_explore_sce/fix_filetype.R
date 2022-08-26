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

lobstr::obj_size(sce)
# 181.85 MB

saveHDF5SummarizedExperiment(sce, here("processed-data", "sce", "sce_DLPFC_annotated"), replace = TRUE)
