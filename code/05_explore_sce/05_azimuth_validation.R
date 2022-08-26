######### Cell Annotation with Azimuth#####
# devtools::install_github("satijalab/seurat-data") # install SeuratData package
# devtools::install_github("satijalab/azimuth", ref = "release/0.4.5") # install Azimuth package

library("SingleCellExperiment")
library("Seurat")
library("Azimuth")
library("SeuratData")
library("patchwork")
library("here")
# library("pheatmap")
# library("scater")
library("sessioninfo")

#### Plot Setup ####
plot_dir <- here("plots", "05_explore_sce", "05_azimuth_validation")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Cell Annotation with RunAzimuth ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

query <- CreateSeuratObject(
    counts = as.matrix(counts(sce)),
    meta.data = data.frame(colData(sce)),
    project = "DLPFC"
)

## Use Azimuth Refrence ##
set.seed(20220818)
query <- RunAzimuth(query, reference = "humancortexref") ## Cell annotation with Azimuth


sce$cellType_azimuth <- query$predicted.subclass
table(query$predicted.subclass)
# Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT      L6 IT L6 IT Car3        L6b
# 6744       9888      21440        109       6775        345       1158       1155        332        875
# Lamp5  Micro-PVM      Oligo        OPC      Pvalb       Sncg        Sst  Sst Chodl        Vip       VLMC
# 802       2507      11384       1836       1513        243       4689         87       4178       1544

## Save
SeuratDisk::SaveH5Seurat(query,
    filename = here("processed-data", "05_explore_sce", "05_azimuth_validation", "sce_DLPFC.h5Seurat"),
    overwrite = TRUE
)

save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('05_azimuth_validation', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript 05_azimuth_validation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R Under development (unstable) (2021-11-06 r81149)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# iocParallel                1.30.3     2022-06-05 [2] Bioconductor
# BiocSingular                1.12.0     2022-04-26 [2] Bioconductor
# bit                         4.0.4      2020-08-04 [2] CRAN (R 4.1.0)
# bit64                       4.0.5      2020-08-30 [2] CRAN (R 4.1.0)
# bitops                      1.0-7      2021-04-24 [2] CRAN (R 4.2.0)
# cellranger                  1.1.0      2016-07-27 [2] CRAN (R 4.1.0)
# cli                         3.3.0      2022-04-25 [2] CRAN (R 4.2.0)
# cluster                     2.1.3      2022-03-28 [3] CRAN (R 4.2.0)
# codetools                   0.2-18     2020-11-04 [3] CRAN (R 4.2.0)
# colorout                  * 1.2-2      2022-04-06 [1] Github (jalvesaq/colorout@79931fd)
# colorspace                  2.0-3      2022-02-21 [2] CRAN (R 4.2.0)
# cowplot                     1.1.1      2020-12-30 [2] CRAN (R 4.2.0)
# crayon                      1.5.1      2022-03-26 [2] CRAN (R 4.2.0)
# curl                        4.3.2      2021-06-23 [2] CRAN (R 4.2.0)
# data.table                  1.14.2     2021-09-27 [2] CRAN (R 4.2.0)
# DBI                         1.1.3      2022-06-18 [2] CRAN (R 4.2.0)
# DelayedArray              * 0.22.0     2022-04-26 [2] Bioconductor
# DelayedMatrixStats          1.18.0     2022-04-26 [2] Bioconductor
# deldir                      1.0-6      2021-10-23 [1] CRAN (R 4.2.0)
# digest                      0.6.29     2021-12-01 [2] CRAN (R 4.2.0)
# dplyr                       1.0.9      2022-04-28 [2] CRAN (R 4.2.0)
# DT                          0.24       2022-08-09 [2] CRAN (R 4.2.0)
# ellipsis                    0.3.2      2021-04-29 [2] CRAN (R 4.2.0)
# fansi                       1.0.3      2022-03-24 [2] CRAN (R 4.2.0)
# fastmap                     1.1.0      2021-01-25 [2] CRAN (R 4.1.0)
# fitdistrplus                1.1-8      2022-03-10 [1] CRAN (R 4.2.0)
# fs                          1.5.2      2021-12-08 [2] CRAN (R 4.2.0)
# future                      1.27.0     2022-07-22 [2] CRAN (R 4.2.0)
# future.apply                1.9.0      2022-04-25 [2] CRAN (R 4.2.0)
# gargle                      1.2.0      2021-07-02 [2] CRAN (R 4.2.0)
# generics                    0.1.3      2022-07-05 [2] CRAN (R 4.2.0)
# GenomeInfoDb              * 1.32.3     2022-08-09 [2] Bioconductor
# GenomeInfoDbData            1.2.8      2022-04-16 [2] Bioconductor
# GenomicRanges             * 1.48.0     2022-04-26 [2] Bioconductor
# ggplot2                     3.3.6      2022-05-03 [2] CRAN (R 4.2.0)
# ggrepel                     0.9.1      2021-01-15 [2] CRAN (R 4.1.0)
# ggridges                    0.5.3      2021-01-08 [2] CRAN (R 4.1.0)
# globals                     0.16.0     2022-08-05 [2] CRAN (R 4.2.0)
# glue                        1.6.2      2022-02-24 [2] CRAN (R 4.2.0)
# goftest                     1.2-3      2021-10-07 [1] CRAN (R 4.2.0)
# googledrive                 2.0.0      2021-07-08 [2] CRAN (R 4.2.0)
# googlesheets4               1.0.1      2022-08-13 [2] CRAN (R 4.2.0)
# gridExtra                   2.3        2017-09-09 [2] CRAN (R 4.1.0)
# gtable                      0.3.0      2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array                 * 1.24.2     2022-08-02 [2] Bioconductor
# hdf5r                       1.3.5      2021-11-15 [2] CRAN (R 4.2.0)
# here                      * 1.0.1      2020-12-13 [1] CRAN (R 4.2.0)
# htmltools                   0.5.3      2022-07-18 [1] CRAN (R 4.2.0)
# htmlwidgets                 1.5.4      2021-09-08 [2] CRAN (R 4.2.0)
# httpuv                      1.6.5      2022-01-05 [2] CRAN (R 4.2.0)
# httr                        1.4.4      2022-08-17 [2] CRAN (R 4.2.0)
# humancortexref.SeuratData * 1.0.0      2022-08-16 [1] local
# ica                         1.0-3      2022-07-08 [1] CRAN (R 4.2.0)
# igraph                      1.3.4      2022-07-19 [2] CRAN (R 4.2.0)
# IRanges                   * 2.30.0     2022-04-26 [2] Bioconductor
# irlba                       2.3.5      2021-12-06 [2] CRAN (R 4.2.0)
# jsonlite                    1.8.0      2022-02-22 [2] CRAN (R 4.2.0)
# KernSmooth                  2.23-20    2021-05-03 [3] CRAN (R 4.2.0)
# later                       1.3.0      2021-08-18 [2] CRAN (R 4.2.0)
# lattice                     0.20-45    2021-09-22 [3] CRAN (R 4.2.0)
# lazyeval                    0.2.2      2019-03-15 [2] CRAN (R 4.1.0)
# leiden                      0.4.2      2022-05-09 [1] CRAN (R 4.2.0)
# lifecycle                   1.0.1      2021-09-24 [2] CRAN (R 4.2.0)
# listenv                     0.8.0      2019-12-05 [2] CRAN (R 4.1.0)
# lmtest                      0.9-40     2022-03-21 [2] CRAN (R 4.2.0)
# magrittr                    2.0.3      2022-03-30 [2] CRAN (R 4.2.0)
# MASS                        7.3-58.1   2022-08-03 [3] CRAN (R 4.2.0)
# Matrix                    * 1.4-1      2022-03-23 [3] CRAN (R 4.2.0)
# MatrixGenerics            * 1.8.1      2022-06-26 [2] Bioconductor
# matrixStats               * 0.62.0     2022-04-19 [2] CRAN (R 4.2.0)
# mgcv                        1.8-40     2022-03-29 [3] CRAN (R 4.2.0)
# mime                        0.12       2021-09-28 [2] CRAN (R 4.2.0)
# miniUI                      0.1.1.1    2018-05-18 [2] CRAN (R 4.1.0)
# munsell                     0.5.0      2018-06-12 [2] CRAN (R 4.1.0)
# nlme                        3.1-159    2022-08-09 [3] CRAN (R 4.2.0)
# parallelly                  1.32.1     2022-07-21 [2] CRAN (R 4.2.0)
# patchwork                 * 1.1.1      2020-12-17 [2] CRAN (R 4.2.0)
# pbapply                     1.5-0      2021-09-16 [1] CRAN (R 4.2.0)
# pillar                      1.8.0      2022-07-18 [1] CRAN (R 4.2.0)
# pkgconfig                   2.0.3      2019-09-22 [2] CRAN (R 4.1.0)
# plotly                      4.10.0     2021-10-09 [2] CRAN (R 4.2.0)
# plyr                        1.8.7      2022-03-24 [2] CRAN (R 4.2.0)
# png                         0.1-7      2013-12-03 [2] CRAN (R 4.1.0)
# polyclip                    1.10-0     2019-03-14 [2] CRAN (R 4.1.0)
# presto                      1.0.0      2022-08-12 [1] Github (immunogenomics/presto@052085d)
# progressr                   0.10.1     2022-06-03 [2] CRAN (R 4.2.0)
# promises                    1.2.0.1    2021-02-11 [2] CRAN (R 4.1.0)
# purrr                       0.3.4      2020-04-17 [2] CRAN (R 4.1.0)
# R6                          2.5.1      2021-08-19 [2] CRAN (R 4.2.0)
# RANN                        2.6.1      2019-01-08 [2] CRAN (R 4.1.0)
# rappdirs                    0.3.3      2021-01-31 [2] CRAN (R 4.1.0)
# RColorBrewer                1.1-3      2022-04-03 [2] CRAN (R 4.2.0)
# Rcpp                        1.0.9      2022-07-08 [2] CRAN (R 4.2.0)
# RcppAnnoy                   0.0.19     2021-07-30 [2] CRAN (R 4.2.0)
# RCurl                       1.98-1.8   2022-07-30 [2] CRAN (R 4.2.0)
# reshape2                    1.4.4      2020-04-09 [2] CRAN (R 4.1.0)
# reticulate                  1.25       2022-05-11 [1] CRAN (R 4.2.0)
# rgeos                       0.5-9      2021-12-15 [1] CRAN (R 4.2.0)
# rhdf5                     * 2.40.0     2022-04-26 [2] Bioconductor
# rhdf5filters                1.8.0      2022-04-26 [2] Bioconductor
# Rhdf5lib                    1.18.2     2022-05-15 [2] Bioconductor
# rlang                       1.0.4      2022-07-12 [2] CRAN (R 4.2.0)
# ROCR                        1.0-11     2020-05-02 [2] CRAN (R 4.1.0)
# rpart                       4.1.16     2022-01-24 [3] CRAN (R 4.2.0)
# rprojroot                   2.0.3      2022-04-02 [2] CRAN (R 4.2.0)
# rsvd                        1.0.5      2021-04-16 [2] CRAN (R 4.2.0)
# Rtsne                       0.16       2022-04-17 [2] CRAN (R 4.2.0)
# S4Vectors                 * 0.34.0     2022-04-26 [2] Bioconductor
# ScaledMatrix                1.4.0      2022-04-26 [2] Bioconductor
# scales                      1.2.0      2022-04-13 [2] CRAN (R 4.2.0)
# scattermore                 0.8        2022-02-14 [1] CRAN (R 4.2.0)
# scry                        1.8.0      2022-04-26 [2] Bioconductor
# sctransform                 0.3.3      2022-01-13 [1] CRAN (R 4.2.0)
# scuttle                     1.6.2      2022-05-15 [2] Bioconductor
# sessioninfo               * 1.2.2      2021-12-06 [2] CRAN (R 4.2.0)
# Seurat                    * 4.1.1      2022-05-02 [1] CRAN (R 4.2.0)
# SeuratData                * 0.2.2      2022-08-12 [1] Github (satijalab/seurat-data@d6a8ce6)
# SeuratDisk                  0.0.0.9020 2022-08-12 [1] Github (mojaveazure/seurat-disk@9b89970)
# SeuratObject              * 4.1.0      2022-05-01 [1] CRAN (R 4.2.0)
# shiny                       1.7.2      2022-07-19 [2] CRAN (R 4.2.0)
# shinyBS                   * 0.61.1     2022-04-17 [1] CRAN (R 4.2.0)
# shinydashboard              0.7.2      2021-09-30 [2] CRAN (R 4.2.0)
# shinyjs                     2.1.0      2021-12-23 [2] CRAN (R 4.2.0)
# SingleCellExperiment      * 1.18.0     2022-04-26 [2] Bioconductor
# sp                        * 1.5-0      2022-06-05 [2] CRAN (R 4.2.0)
# sparseMatrixStats           1.8.0      2022-04-26 [2] Bioconductor
# spatstat.core               2.4-4      2022-05-18 [1] CRAN (R 4.2.0)
# spatstat.data               2.2-0      2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.geom               2.4-0      2022-03-29 [1] CRAN (R 4.2.0)
# spatstat.random             2.2-0      2022-03-30 [1] CRAN (R 4.2.0)
# spatstat.sparse             2.1-1      2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.utils              2.3-1      2022-05-06 [1] CRAN (R 4.2.0)
# stringi                     1.7.8      2022-07-11 [2] CRAN (R 4.2.0)
# stringr                     1.4.0      2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment      * 1.26.1     2022-04-29 [2] Bioconductor
# survival                    3.4-0      2022-08-09 [3] CRAN (R 4.2.0)
# tensor                      1.5        2012-05-05 [1] CRAN (R 4.2.0)
# tibble                      3.1.8      2022-07-22 [2] CRAN (R 4.2.0)
# tidyr                       1.2.0      2022-02-01 [2] CRAN (R 4.2.0)
# tidyselect                  1.1.2      2022-02-21 [2] CRAN (R 4.2.0)
# utf8                        1.2.2      2021-07-24 [2] CRAN (R 4.2.0)
# uwot                        0.1.13     2022-08-16 [2] CRAN (R 4.2.0)
# vctrs                       0.4.1      2022-04-13 [2] CRAN (R 4.2.0)
# viridisLite                 0.4.0      2021-04-13 [2] CRAN (R 4.2.0)
# withr                       2.5.0      2022-03-03 [2] CRAN (R 4.2.0)
# xtable                      1.8-4      2019-04-21 [2] CRAN (R 4.1.0)
# XVector                     0.36.0     2022-04-26 [2] Bioconductor
# zlibbioc                    1.42.0     2022-04-26 [2] Bioconductor
# zoo                         1.8-10     2022-04-15 [2] CRAN (R 4.2.0)
#
# [1] /users/lhuuki/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#
#
