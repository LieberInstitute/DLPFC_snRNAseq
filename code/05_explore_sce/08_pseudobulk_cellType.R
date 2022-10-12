
library("SingleCellExperiment")
library("scuttle")
library("edgeR")
library("here")
library("sessioninfo")

args <- commandArgs(trailingOnly = TRUE)
cellType <- args[1]
message("\n#### Pseudobulking: ", dataset, " ####")

## load data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

data_dir <- here("processed-data", "05_explore_sce", "08_pseudobulk_cellTypes")
if(!dir.exists(data_dir)) dir.create(data_dir)

#### Pseudobulk cellType - hc ####
## remove an NAs
sce <- sce[, !is.na(sce$cellType_layer)]

message(Sys.time(), "- Pseudobulking")

sce_pseudo <- aggregateAcrossCells(sce,
                                 ids = paste0(sce$Sample ,"-", sce[[cellType]]))

## Any very small groups? (n < 10)
message("Any small groups?")
table(sce_pseudo$ncells < 10)

sce_pseudo <- sce_pseudo[, sce_pseudo$ncells >= 10]

## Need to normalize
message(Sys.time(), "- Normalize")

logcounts(sce_pseudo) <-
  edgeR::cpm(edgeR::calcNormFactors(sce_pseudo),
             log = TRUE,
             prior.count = 1
  )

message(Sys.time(), "- Save Data")
save(sce_pseudo, file = here(data_dir, file = paste0("sce_pseudo-", cellType, ".Rdata")))

# sgejobs::job_single('08_pseudobulk_cellType_layer', create_shell = TRUE, memory = '25G', command = "Rscript 08_pseudobulk_cellType.R cellType_hc")
# sgejobs::job_single('08_pseudobulk_cellType_hc', create_shell = TRUE, memory = '25G', command = "Rscript 08_pseudobulk_cellType.R cellType_layer")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
