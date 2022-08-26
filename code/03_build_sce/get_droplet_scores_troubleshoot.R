## Based on
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library("SingleCellExperiment")
library("DropletUtils")
# library("BiocParallel")
library("scuttle")
library("tidyverse")
library("here")
library("sessioninfo")

## get user supplied sample + lower cutoff
args <- commandArgs(trailingOnly = TRUE)
user_sample <- args[[1]]
user_lower <- as.integer(args[[2]])

#### Load & Subset raw data ####
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)

stopifnot(user_sample %in% sce$Sample)

message("Running Sample: ", user_sample)

sce <- sce[, sce$Sample == user_sample]
message("ncol:", ncol(sce))

#### Run barcodeRanks to find knee ####

bcRanks <- barcodeRanks(sce, fit.bounds = c(10, 1e3))

knee_lower <- metadata(bcRanks)$knee + 100
message(
    "'Second knee point' = ", metadata(bcRanks)$knee, "\n",
    "knee_lower = ", knee_lower, "\n",
    "user defined lower = ", user_lower
)

#### Run emptyDrops w/ knee + 100 ####
set.seed(100)
message("Starting emptyDrops")
Sys.time()
e.out <- DropletUtils::emptyDrops(
    sce,
    niters = 30000,
    lower = user_lower
    # ,
    # BPPARAM = BiocParallel::MulticoreParam(4)
)
message("Done - saving data")
Sys.time()

save(e.out, file = here("processed-data", "03_build_sce", "droplet_scores_troubleshoot", paste0("droplet_scores_", user_sample, ".Rdata")))

#### QC Plots ####
message("QC check")
FDR_cutoff <- 0.001
addmargins(table(Signif = e.out$FDR <= FDR_cutoff, Limited = e.out$Limited, useNA = "ifany"))

n_cell_anno <- paste("Non-empty:", sum(e.out$FDR < FDR_cutoff, na.rm = TRUE))
message(n_cell_anno)

my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

droplet_elbow_plot <- as.data.frame(bcRanks) %>%
    add_column(FDR = e.out$FDR) %>%
    ggplot(aes(x = rank, y = total, color = FDR < FDR_cutoff)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = metadata(bcRanks)$knee, linetype = "dotted", color = "gray") +
    annotate("text", x = 10, y = metadata(bcRanks)$knee, label = "Second Knee", vjust = -1, color = "gray") +
    geom_hline(yintercept = knee_lower, linetype = "dashed") +
    annotate("text", x = 10, y = knee_lower, label = "Knee est 'lower'", vjust = -0.5) +
    geom_hline(yintercept = user_lower, linetype = "dashed", color = "red") +
    annotate("text", x = 10, y = user_lower, label = "User defined 'lower'", vjust = -0.5, color = "red") +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    labs(
        x = "Barcode Rank",
        y = "Total UMIs",
        title = paste("Sample", user_sample),
        subtitle = n_cell_anno,
        color = paste("FDR <", FDR_cutoff)
    ) +
    my_theme +
    theme(legend.position = "bottom")

# droplet_scatter_plot <- as.data.frame(e) %>%
#   ggplot(aes(x = Total, y = -LogProb, color = FDR < FDR_cutoff)) +
#   geom_point(alpha = 0.5, size = 1) +
#   labs(x = "Total UMIs", y = "-Log Probability",
#        color = paste("FDR <", FDR_cutoff)) +
#   my_theme+
#   theme(legend.position = "bottom")
# # print(droplet_elbow_plot/droplet_scatter_plot)
# ggsave(droplet_elbow_plot/droplet_scatter_plot, filename = here("plots","03_build_sce", "droplet_qc_png",paste0("droplet_qc_",sample,".png")))

ggsave(droplet_elbow_plot, filename = here("plots", "03_build_sce", "droplet_qc", "droplet_knee_plots", paste0("droplet_qc_", user_sample, ".png")))


# sgejobs::job_single('get_droplet_scores_troubleshoot', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_droplet_scores_troubleshoot.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# [1] "Reproducibility information:"
# [1] "2022-05-05 15:37:47 EDT"
# user   system  elapsed
# 1487.874   20.197 1510.183
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
# date     2022-05-05
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.1.2)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocParallel           1.28.3   2021-12-09 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# broom                  0.8.0    2022-04-13 [2] CRAN (R 4.1.2)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                * 1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.1.2)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.5.0    2022-04-15 [2] CRAN (R 4.1.2)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.1.2)
# httr                   1.4.3    2022-05-04 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.1.2)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [1] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
# ragg                   1.2.2    2022-02-21 [2] CRAN (R 4.1.2)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.1.2)
# readxl                 1.4.0    2022-03-28 [2] CRAN (R 4.1.2)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.1.1)
# rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
# scales                 1.2.0    2022-04-13 [2] CRAN (R 4.1.2)
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.1.2)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# systemfonts            1.0.4    2022-02-11 [2] CRAN (R 4.1.2)
# textshaping            0.3.6    2021-10-13 [2] CRAN (R 4.1.2)
# tibble               * 3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.1.0)
# tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/lhuuki/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
