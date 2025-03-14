## Louise Huuki-Myers, MArch 2025
## Run spatialLIBD::registration_wrapper to get cell type modeling for Layer annotations

## Required libraries
library("here")
library("sessioninfo")
library("SingleCellExperiment")
library("spatialLIBD")
library("HDF5Array")

# Import command-line parameters

data_dir <- here("processed-data", "05_explore_sce", "10_model_cellType_layer")
if(!dir.exists(data_dir)) dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

#### Load the data ####
message(Sys.time(), " - Load Data")
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## make APOE syntatic

modeling_results <-registration_wrapper(
  sce = sce,
  var_registration = "cellType_layer",
  var_sample_id = "Sample",
  covars = c("Position", "age", "sex"),
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

message(Sys.time(), " - Saving Data")
saveRDS(modeling_results, file = here(data_dir, "DLPFCsn_modeling_results-cellType_layer.rds"))

# slurmjobs::job_single('10_model_cellType_layer', create_shell = TRUE, memory = '100G', command = "Rscript 10_model_cellType_layer.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


