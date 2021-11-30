# ssh j
# screen -S snRNAseq
# qrsh -l bluejay,mem_free=25G,h_vmem=25G -pe local 4

## Based on
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library("SingleCellExperiment")
library("DropletUtils")
library("BiocParallel")
library("scuttle")
library("here")
library("sessioninfo")

## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)

## Find empty droplets
Sys.time()
e.out <- DropletUtils::emptyDrops(
    sce,
    niters = 30000,
    BPPARAM = BiocParallel::MulticoreParam(4)
)
Sys.time()
# [1] "2021-11-30 15:48:22 EST"
# [1] "2021-11-30 16:43:00 EST"

## Check empty droplet results
addmargins(table(Signif = e.out$FDR <= 0.001, Limited = e.out$Limited, useNA = "ifany"))
#        Limited
# Signif     FALSE     TRUE     <NA>      Sum
#   FALSE  1275410        0        0  1275410
#   TRUE     33021   136704        0   169725
#   <NA>         0        0 25160671 25160671
#   Sum    1308431   136704 25160671 26605806
  
## Eliminate empty droplets
sce <- sce[, which(e.out$FDR <= 0.001)]

## Compute QC metrics
sce <- addPerCellQC(
    sce,
    subsets = list(Mito = which(seqnames(sce) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)

## Save for later
save(sce, file = here::here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info  ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  hash: flag: Rwanda, pick, bank
#
#  setting  value
#  version  R version 4.1.2 Patched (2021-11-04 r81138)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-11-30
#  pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date (UTC) lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  beachmat               2.10.0   2021-10-26 [2] Bioconductor
#  Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
#  BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
#  BiocParallel         * 1.28.2   2021-11-25 [2] Bioconductor
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
#  colorout               1.2-2    2021-11-02 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
#  DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
#  digest                 0.6.28   2021-09-23 [2] CRAN (R 4.1.2)
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  dqrng                  0.3.0    2021-05-01 [1] CRAN (R 4.1.2)
#  DropletUtils         * 1.14.1   2021-11-08 [1] Bioconductor
#  edgeR                  3.36.0   2021-10-26 [2] Bioconductor
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
#  GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
#  GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
#  GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
#  ggplot2                3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  glue                   1.5.1    2021-11-30 [2] CRAN (R 4.1.2)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.2)
#  htmltools              0.5.2    2021-08-25 [2] CRAN (R 4.1.2)
#  htmlwidgets            1.5.4    2021-09-08 [2] CRAN (R 4.1.2)
#  httpuv                 1.6.3    2021-09-09 [2] CRAN (R 4.1.2)
#  IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  later                  1.3.0    2021-08-18 [2] CRAN (R 4.1.2)
#  lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
#  lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
#  limma                  3.50.0   2021-10-26 [2] Bioconductor
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.2)
#  MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
#  matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                 1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
#  R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
#  Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
#  rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
#  rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
#  Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
#  rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
#  rmote                  0.3.4    2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scuttle              * 1.4.0    2021-10-26 [1] Bioconductor
#  servr                  0.23     2021-08-11 [1] CRAN (R 4.1.2)
#  sessioninfo          * 1.2.1    2021-11-02 [2] CRAN (R 4.1.2)
#  SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
#  sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
#  SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
#  tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  xfun                   0.28     2021-11-04 [2] CRAN (R 4.1.2)
#  XVector                0.34.0   2021-10-26 [2] Bioconductor
#  zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
