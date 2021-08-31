## Locate the CellRanger metrics from Tran et al

library("here")
library("sessioninfo")

## Load common functions
source(here("code", "02_cellranger_metrics", "00_metrics_functions.R"))

## Locate the metrics summary files first
metrics_files <-
    unlist(sapply(
        c(
            "/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/",
            "/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/Feb2021"
        ),
        find_metrics_csv
    ), use.names = FALSE)

## Read in the metrics
tran_metrics_raw <- lapply(metrics_files, read.csv)
table(sapply(tran_metrics_raw, ncol))
# 19 20
# 30  6

metrics_files[sapply(tran_metrics_raw, ncol) == 20]
tran_metrics_raw[sapply(tran_metrics_raw, ncol) == 20]

tran_metrics_raw <- mapply(function(x, y) {
    cbind(x, metrics_csv = y)
}, tran_metrics_raw, metrics_files)

## Merge
tran_metrics <- Reduce(function(...) merge(..., all = TRUE), tran_metrics_raw)
dim(tran_metrics)
# [1] 36 21

## Simplify to numbers
tran_metrics <- cbind(
    metrics_to_numbers(tran_metrics[, -which(colnames(tran_metrics) == "metrics_csv")]),
    metrics_csv = tran_metrics$metrics_csv
)

## Add information about the samples
tran_metrics$set <- ifelse(grepl("Feb2021", tran_metrics$metrics_csv), "Tran2020", "Tran2021")

## Save for later
dir.create(here("processed-data", "02_cellranger_metrics"), showWarnings = FALSE)
save(tran_metrics, file = here("processed-data", "02_cellranger_metrics", "tran_metrics.Rdata"))
write.csv(tran_metrics, file = here("processed-data", "02_cellranger_metrics", "tran_metrics.csv"))
