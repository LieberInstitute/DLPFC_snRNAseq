## Based on
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
# library("jaffelab")
library("here")
library("sessioninfo")

## get sample #
args = commandArgs(trailingOnly=TRUE)
sample_i = as.integer(args[[1]])

message("sample_i :", sample_i)
## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)

samples <- unique(sce$Sample)
sample_run <- samples[[sample_i]]
message("Running Sample: ", sample_run, " (", sample_i, "/", length(samples),")")

sce <- sce[, sce$Sample == sample_run]
message("ncol:", ncol(sce))

set.seed(100)
message("Starting emptyDrops")
Sys.time()
e.out <- DropletUtils::emptyDrops(
  sce,
  niters = 30000
  # ,
  # BPPARAM = BiocParallel::MulticoreParam(4)
)
message("Done - saving data")
Sys.time()

save(e.out, file = here("processed-data", "03_build_sce", "droplet_scores",paste0("droplet_scores_", sample_run,".Rdata")))

# sgejobs::job_single('get_droplet_scores', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_droplet_scores.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
