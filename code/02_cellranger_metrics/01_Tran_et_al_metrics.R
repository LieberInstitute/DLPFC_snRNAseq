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
dir.create(here("processed-data", "02_cellranger_metrics"),
    showWarnings = FALSE)
save(
    tran_metrics,
    file = here(
        "processed-data",
        "02_cellranger_metrics",
        "tran_metrics.Rdata"
    )
)
write.csv(
    tran_metrics,
    file = here(
        "processed-data",
        "02_cellranger_metrics",
        "tran_metrics.csv"
    )
)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.0 Patched (2021-05-18 r80330)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-08-31
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version date       lib source
#  assertthat    0.2.1   2019-03-21 [2] CRAN (R 4.1.0)
#  cli           3.0.1   2021-07-17 [2] CRAN (R 4.1.0)
#  colorout      1.2-2   2021-05-25 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace    2.0-2   2021-06-24 [2] CRAN (R 4.1.0)
#  crayon        1.4.1   2021-02-08 [2] CRAN (R 4.1.0)
#  DBI           1.1.1   2021-01-15 [2] CRAN (R 4.1.0)
#  digest        0.6.27  2020-10-24 [2] CRAN (R 4.1.0)
#  dplyr         1.0.7   2021-06-18 [2] CRAN (R 4.1.0)
#  ellipsis      0.3.2   2021-04-29 [2] CRAN (R 4.1.0)
#  fansi         0.5.0   2021-05-25 [2] CRAN (R 4.1.0)
#  generics      0.1.0   2020-10-31 [2] CRAN (R 4.1.0)
#  ggplot2       3.3.5   2021-06-25 [2] CRAN (R 4.1.0)
#  glue          1.4.2   2020-08-27 [2] CRAN (R 4.1.0)
#  gtable        0.3.0   2019-03-25 [2] CRAN (R 4.1.0)
#  here        * 1.0.1   2020-12-13 [1] CRAN (R 4.1.0)
#  htmltools     0.5.1.1 2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets   1.5.3   2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv        1.6.1   2021-05-07 [2] CRAN (R 4.1.0)
#  jsonlite      1.7.2   2020-12-09 [2] CRAN (R 4.1.0)
#  later         1.2.0   2021-04-23 [2] CRAN (R 4.1.0)
#  lattice       0.20-44 2021-05-02 [3] CRAN (R 4.1.0)
#  lifecycle     1.0.0   2021-02-15 [2] CRAN (R 4.1.0)
#  magrittr      2.0.1   2020-11-17 [2] CRAN (R 4.1.0)
#  munsell       0.5.0   2018-06-12 [2] CRAN (R 4.1.0)
#  pillar        1.6.2   2021-07-29 [2] CRAN (R 4.1.0)
#  pkgconfig     2.0.3   2019-09-22 [2] CRAN (R 4.1.0)
#  png           0.1-7   2013-12-03 [2] CRAN (R 4.1.0)
#  promises      1.2.0.1 2021-02-11 [2] CRAN (R 4.1.0)
#  purrr         0.3.4   2020-04-17 [2] CRAN (R 4.1.0)
#  R6            2.5.0   2020-10-28 [2] CRAN (R 4.1.0)
#  Rcpp          1.0.7   2021-07-07 [2] CRAN (R 4.1.0)
#  rlang         0.4.11  2021-04-30 [2] CRAN (R 4.1.0)
#  rmote         0.3.4   2021-05-25 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot     2.0.2   2020-11-15 [2] CRAN (R 4.1.0)
#  scales        1.1.1   2020-05-11 [2] CRAN (R 4.1.0)
#  servr         0.22    2021-04-14 [1] CRAN (R 4.1.0)
#  sessioninfo * 1.1.1   2018-11-05 [2] CRAN (R 4.1.0)
#  tibble        3.1.3   2021-07-23 [2] CRAN (R 4.1.0)
#  tidyselect    1.1.1   2021-04-30 [2] CRAN (R 4.1.0)
#  utf8          1.2.2   2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs         0.3.8   2021-04-29 [2] CRAN (R 4.1.0)
#  withr         2.4.2   2021-04-18 [2] CRAN (R 4.1.0)
#  xfun          0.25    2021-08-06 [2] CRAN (R 4.1.0)
#
# [1] /users/lcollado/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
