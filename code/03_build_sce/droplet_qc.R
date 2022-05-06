library("SingleCellExperiment")
library("scuttle")
library("tidyverse")
library("patchwork")
library("here")
library("sessioninfo")

## plot set up
my_theme <- theme_bw() +
  theme(text = element_text(size=15))


## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)
length(table(sce$Sample))

droplet_score_fn <- list.files(here("processed-data", "03_build_sce","droplet_scores"),
                               full.names = TRUE)

names(droplet_score_fn)  <- gsub("droplet_scores_|.Rdata","",basename(droplet_score_fn))

e.out <- lapply(droplet_score_fn, function(x) get(load(x)))

## check out n empty with boxplot
FDR_cutoff <- 0.001

drop_summary <- stack(map_int(e.out, nrow)) %>%
  rename(total_n = values) %>% 
  left_join(stack(map_int(e.out, ~sum(.x$FDR < FDR_cutoff, na.rm = TRUE))) %>%
              rename(non_empty = values)) %>%
  select(Sample = ind, total_n, non_empty)

write_csv(drop_summary, file = here("processed-data", "03_build_sce","drop_summary.csv"))

drop_summary %>%
  arrange(non_empty)

summary(drop_summary$non_empty)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2989    3682    4134    4461    5270    6269 

## compare to default drop empty
drop_summary_defualt <- read.csv(here("processed-data", "03_build_sce","drop_summary_default.csv"))

drop_summary %>% 
  left_join(drop_summary_defualt) %>% 
  mutate(d = non_empty_default - non_empty) %>%
  arrange(d)

#         Sample total_n non_empty non_empty_default     d
# 1   Br2743_ant  655424      3198              2975  -223
# 2   Br6432_ant  736110      3432              3402   -30
# 3   Br6471_ant 1148322      3622              3594   -28
# 4  Br8492_post  965917      2989              3678   689
# 5   Br2720_mid 1152093      3743              4435   692
# 6   Br6522_mid  604278      4035              4897   862
# 7   Br2743_mid 1495673      3426              4837  1411
# 8  Br6522_post  681586      4295              6066  1771
# 9   Br3942_ant 1692299      6269             52946 46677
# 10  Br6471_mid 1653348      5352             56272 50920
# 11  Br8492_mid 1302254      5189             56478 51289
# 12 Br2720_post  922885      5950             57625 51675
# 13  Br8325_ant 1855444      5951             58687 52736
# 14  Br8667_ant 1875106      5806             70117 64311
# 15  Br3942_mid 2209670      4309             74365 70056
# 16  Br8325_mid 2192531      4990             77615 72625
# 17  Br6423_ant 1521781      3999             88179 84180
# 18  Br8667_mid 2154338      4134             88660 84526
# 19 Br6423_post 1786747      4067             91276 87209

drop_barplot <- drop_summary %>%
  mutate(empty = total_n - non_empty) %>%
  select(-total_n) %>%
  pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
  ggplot(aes(x = Sample, y = n_drop, fill = drop_type))+
  geom_col() +
  scale_y_continuous(trans='log10') +
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust = 1))

ggsave(drop_barplot, filename = here("plots","03_build_sce", "drop_barplot.png"), width = 9)

## Check empty droplet results
map(e.out, ~addmargins(table(Signif = .x$FDR <= FDR_cutoff, Limited = .x$Limited, useNA = "ifany")))

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce),]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1] 36601 84756

## Compute QC metrics
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = which(seqnames(sce) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

## Save for later
save(sce, file = here::here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

# sgejobs::job_single('droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript droplet_qc.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R Under development (unstable) (2021-11-06 r81149)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-05-06
# pandoc   2.11.0.4 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/bin/pandoc
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# backports              1.4.1    2021-12-13 [2] CRAN (R 4.2.0)
# beachmat               2.12.0   2022-04-26 [2] Bioconductor
# Biobase              * 2.56.0   2022-04-26 [2] Bioconductor
# BiocGenerics         * 0.42.0   2022-04-26 [2] Bioconductor
# BiocParallel           1.30.0   2022-04-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.2.0)
# broom                  0.8.0    2022-04-13 [2] CRAN (R 4.2.0)
# cellranger             1.1.0    2016-07-27 [2] CRAN (R 4.1.0)
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.2.0)
# colorout             * 1.2-2    2022-04-06 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.2.0)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.2.0)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.2.0)
# dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
# DelayedArray           0.22.0   2022-04-26 [2] Bioconductor
# DelayedMatrixStats     1.18.0   2022-04-26 [2] Bioconductor
# dplyr                * 1.0.9    2022-04-28 [2] CRAN (R 4.2.0)
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.2.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.2.0)
# forcats              * 0.5.1    2021-01-27 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.2.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.2.0)
# GenomeInfoDb         * 1.32.1   2022-04-28 [2] Bioconductor
# GenomeInfoDbData       1.2.8    2022-04-16 [2] Bioconductor
# GenomicRanges        * 1.48.0   2022-04-26 [2] Bioconductor
# ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.2.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.2.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# haven                  2.5.0    2022-04-15 [2] CRAN (R 4.2.0)
# here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.2.0)
# hms                    1.1.1    2021-09-26 [2] CRAN (R 4.2.0)
# httr                   1.4.3    2022-05-04 [2] CRAN (R 4.2.0)
# IRanges              * 2.30.0   2022-04-26 [2] Bioconductor
# jsonlite               1.8.0    2022-02-22 [2] CRAN (R 4.2.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.2.0)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.2.0)
# lubridate              1.8.0    2021-10-07 [2] CRAN (R 4.2.0)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.2.0)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.2.0)
# MatrixGenerics       * 1.8.0    2022-04-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.2.0)
# modelr                 0.1.8    2020-05-19 [2] CRAN (R 4.1.0)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# patchwork            * 1.1.1    2020-12-17 [2] CRAN (R 4.2.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.2.0)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                * 0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.2.0)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.2.0)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.2.0)
# readr                * 2.1.2    2022-01-30 [2] CRAN (R 4.2.0)
# readxl                 1.4.0    2022-03-28 [2] CRAN (R 4.2.0)
# reprex                 2.0.1    2021-08-05 [2] CRAN (R 4.2.0)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.2.0)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.2.0)
# rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
# rvest                  1.0.2    2021-10-16 [2] CRAN (R 4.2.0)
# S4Vectors            * 0.34.0   2022-04-26 [2] Bioconductor
# scales                 1.2.0    2022-04-13 [2] CRAN (R 4.2.0)
# scuttle              * 1.6.0    2022-04-26 [2] Bioconductor
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.2.0)
# SingleCellExperiment * 1.18.0   2022-04-26 [2] Bioconductor
# sparseMatrixStats      1.8.0    2022-04-26 [2] Bioconductor
# stringi                1.7.6    2021-11-29 [2] CRAN (R 4.2.0)
# stringr              * 1.4.0    2019-02-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.26.1   2022-04-29 [2] Bioconductor
# tibble               * 3.1.7    2022-05-03 [2] CRAN (R 4.2.0)
# tidyr                * 1.2.0    2022-02-01 [2] CRAN (R 4.2.0)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.2.0)
# tidyverse            * 1.3.1    2021-04-15 [2] CRAN (R 4.2.0)
# tzdb                   0.3.0    2022-03-28 [2] CRAN (R 4.2.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.2.0)
# vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.2.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.2.0)
# xml2                   1.3.3    2021-11-30 [2] CRAN (R 4.2.0)
# XVector                0.36.0   2022-04-26 [2] Bioconductor
# zlibbioc               1.42.0   2022-04-26 [2] Bioconductor
# 
# [1] /users/lhuuki/R/devel
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-devel/R/devel/lib64/R/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
