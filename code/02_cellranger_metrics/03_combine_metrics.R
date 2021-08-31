## Combine our metrics with those from Tran et al

library("here")
library("sessioninfo")

## Load previous files
load(here(
    "processed-data",
    "02_cellranger_metrics",
    "tran_metrics.Rdata"
),
    verbose = TRUE)

load(
    here(
        "processed-data",
        "02_cellranger_metrics",
        "cellranger_metrics.Rdata"
    ),
    verbose = TRUE
)

## Combine metrics
all_metrics <- merge(tran_metrics, cellranger_metrics, all = TRUE)
dim(all_metrics)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
