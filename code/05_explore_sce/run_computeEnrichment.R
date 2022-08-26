library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")
# library("rafalib")
# library("scuttle")
# library(limma)
# library("RColorBrewer")
# library(lattice)
# library(edgeR)
library("here")
library("sessioninfo")

## Load computeEnrichment
source("utils.R")

## get sample i
args <- commandArgs(trailingOnly = TRUE)
var_oi_input <- args[[1]]

# load sc data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## Edit colData
rownames(sce) <- rowData(sce)$gene_id # have to make row names of object the ensembl id instead of gene names
colData(sce)$Position <- as.factor(colData(sce)$Position)
colData(sce)$age <- as.numeric(colData(sce)$age)
colData(sce)$sex <- as.factor(colData(sce)$sex)
colnames(colData(sce))[1] <- "sample_id"

stopifnot(var_oi_input %in% colnames(colData(sce)))

message("Running Enrichment for: ", var_oi_input)

results_specificity <- computeEnrichment(sce, var_oi = var_oi_input, covars = c("Position", "age", "sex"))


message("Done! saving - ", Sys.time())
save(results_specificity, file = here("processed-data", "05_explore_sce", "enrichment", paste0("enrichment_", var_oi_input, ".Rdata")))


# sgejobs::job_single('run_computeEnrichment_hc', create_shell = TRUE, memory = '50G', command = "Rscript run_computeEnrichment.R cellType_hc")
# sgejobs::job_single('run_computeEnrichment_azimuth', create_shell = TRUE, memory = '50G', command = "Rscript run_computeEnrichment.R cellType_azimuth")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
