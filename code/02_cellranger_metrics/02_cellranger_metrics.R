## Locate the CellRanger metrics

library("here")
library("sessioninfo")

## Load common functions
source(here("code", "02_cellranger_metrics", "00_metrics_functions.R"))

## Locate the metrics summary files first
metrics_files <- find_metrics_csv(here("processed-data", "cellranger"))

## Read in the metrics
cellranger_metrics_raw <- lapply(metrics_files, read.csv)
table(sapply(cellranger_metrics_raw, ncol))
# 19
#  3


## Merge
cellranger_metrics <- do.call(rbind, cellranger_metrics_raw)
dim(cellranger_metrics)
# [1]  3 19

## Simplify to numbers
cellranger_metrics <- metrics_to_numbers(cellranger_metrics)

## Add information about the samples
tmp <- read.table(here("raw-data", "sample_libs_round0-5.tsv"), header = FALSE, row.names = 1)

sample_info <- data.frame(
  row.names = rownames(tmp),
  region_short = tolower(gsub("DLPFC_", "", tmp$V2)),
  subject = tmp$V3,
  round = tmp$V5
)
## match tables
sample_info <- sample_info[rownames(cellranger_metrics),]

cellranger_metrics$metrics_csv <- metrics_files
cellranger_metrics$set <- sample_info$round

## Save for later
dir.create(here("processed-data", "02_cellranger_metrics"),
    showWarnings = FALSE
)
save(
    cellranger_metrics,
    file = here(
        "processed-data",
        "02_cellranger_metrics",
        "cellranger_metrics.Rdata"
    )
)
write.csv(
    cellranger_metrics,
    file = here(
        "processed-data",
        "02_cellranger_metrics",
        "cellranger_metrics.csv"
    )
)

# sgejobs::job_single('cellranger_metrics', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 02_cellranger_metrics.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## Same as 01_Tran_et_al_metrics.R
