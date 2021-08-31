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
cellranger_metrics$metrics_csv <- metrics_files
cellranger_metrics$set <- ifelse(grepl("1c-k|2c-k|3c-k", cellranger_metrics$metrics_csv), "round1", NA)


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

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

## Same as 01_Tran_et_al_metrics.R
