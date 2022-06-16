library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("HDF5Array")
library("here")
library("sessioninfo")

## Load HD5F
sce <- loadHDF5SummarizedExperiment(dir=here("processed-data", "03_build_sce", "hdf5_sce"), prefix="")
dim(sce)

#### DEviance Feat Selection ####
set.seed(606)

message("running Deviance Feat. Selection - ", Sys.time())
sce <- devianceFeatureSelection(sce,
                                assay="counts", fam="binomial", sorted=F,
                                batch=as.factor(sce$round))

# This temp file just used for getting batch-corrected components (drops a variety of entries)

pdf(here("plots","03_build_sce","normalize1", "binomial_deviance.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
dev.off()

message("running nullResiduals - ", Sys.time())
sce <- nullResiduals(sce, assay="counts", fam="binomial",  # default params
                     type="deviance")


hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing=T)][1:2000]

message("running PCA - ", Sys.time())
sce_uncorrected <- runPCA(sce, exprs_values="binomial_deviance_residuals",
                      subset_row=hdgs.hb, ncomponents=100,
                      name="GLMPCA_approx",
                      BSPARAM=BiocSingular::IrlbaParam())

## Why multi match as well?
# sce <- multiBatchNorm(sce, batch=sce$sample_short)

#### Reduce Dims####

message("running TSNE - ", Sys.time())
sce_uncorrected <- runTSNE(sce_uncorrected, dimred="GLMPCA_approx")

message("running UMAP - ", Sys.time())
sce_uncorrected <- runUMAP(sce_uncorrected, dimred="GLMPCA_approx")

message("Saving Data - ", Sys.time())
saveHDF5SummarizedExperiment(sce_uncorrected, dir = here("processed-data", "03_build_sce","sce_uncorrected"))
# sgejobs::job_single('normalize_step1_glm', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript normalize_step1_glm.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2022-06-10 00:00:58 EDT"
# user    system   elapsed 
# 10407.244   397.344 10828.361 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-06-10
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster                1.4.0    2021-10-26 [2] Bioconductor
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
# cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
# codetools              0.2-18   2020-11-04 [3] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray         * 0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# dplyr                  1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.1.2)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array            * 1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.3.1    2022-04-20 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [1] CRAN (R 4.1.2)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# Matrix               * 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RcppAnnoy              0.0.19   2021-07-30 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.7 2022-06-09 [2] CRAN (R 4.1.2)
# rhdf5                * 2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
# RSpectra               0.16-1   2022-04-24 [2] CRAN (R 4.1.2)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# Rtsne                  0.16     2022-04-17 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.2.0    2022-04-13 [2] CRAN (R 4.1.2)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scran                * 1.22.1   2021-11-14 [2] Bioconductor
# scry                 * 1.6.0    2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# uwot                   0.1.11   2021-12-02 [2] CRAN (R 4.1.2)
# vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
