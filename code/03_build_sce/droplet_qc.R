library("SingleCellExperiment")
library("scuttle")
library("scran")
library("scater")
library("scDblFinder")
library("jaffelab")
library("batchelor")
library("tidyverse")
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

#### Compute QC metrics ####
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = which(seqnames(sce) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)


#### Check for low quality nuc ####
## High mito
# sce$high.mito.sample ## standard name?
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher", batch = sce$Sample)
table(sce$high_mito)
# FALSE  TRUE 
# 78197  6559 
table(sce$high_mito, sce$Sample)


## low library size
# sce$qc.lib ## standard name?
sce$low_sum <- isOutlier(sce$sum, log=TRUE, type="lower", batch = sce$Sample)
table(sce$low_sum)
# FALSE  TRUE 
# 84294   462


## low detected features
# sce$qc.detected
sce$low_detected <- isOutlier(sce$detected, log=TRUE, type="lower", batch = sce$Sample)
table(sce$low_detected)
# FALSE  TRUE 
# 83513  1243

## All low sum are also low detected
table(sce$low_sum, sce$low_detected)
#       FALSE  TRUE
# FALSE 83513   781
# TRUE      0   462

## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected

table(sce$discard_auto)
# FALSE  TRUE 
# 77604  7152

## discard 8% of nuc
100*sum(sce$discard_auto)/ncol(sce)
# [1] 8.438341

(qc_t <- addmargins(table(sce$Sample, sce$discard_auto)))
#             FALSE  TRUE   Sum
# Br2720_mid   3101   642  3743
# Br2720_post  5911    39  5950
# Br2743_ant   2861   337  3198
# Br2743_mid   2723   703  3426
# Br3942_ant   5205  1064  6269
# Br3942_mid   4282    27  4309
# Br6423_ant   3898   101  3999
# Br6423_post  4067     0  4067
# Br6432_ant   3059   373  3432
# Br6471_ant   3212   410  3622
# Br6471_mid   4724   628  5352
# Br6522_mid   4004    31  4035
# Br6522_post  4275    20  4295
# Br8325_ant   4707  1244  5951
# Br8325_mid   4020   970  4990
# Br8492_mid   4997   192  5189
# Br8492_post  2661   328  2989
# Br8667_ant   5774    32  5806
# Br8667_mid   4123    11  4134
# Sum         77604  7152 84756

round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#             FALSE  TRUE   Sum
# Br2720_mid   82.8  17.2 100.0
# Br2720_post  99.3   0.7 100.0
# Br2743_ant   89.5  10.5 100.0
# Br2743_mid   79.5  20.5 100.0 *
# Br3942_ant   83.0  17.0 100.0
# Br3942_mid   99.4   0.6 100.0
# Br6423_ant   97.5   2.5 100.0
# Br6423_post 100.0   0.0 100.0
# Br6432_ant   89.1  10.9 100.0
# Br6471_ant   88.7  11.3 100.0
# Br6471_mid   88.3  11.7 100.0
# Br6522_mid   99.2   0.8 100.0
# Br6522_post  99.5   0.5 100.0
# Br8325_ant   79.1  20.9 100.0 *
# Br8325_mid   80.6  19.4 100.0 *
# Br8492_mid   96.3   3.7 100.0
# Br8492_post  89.0  11.0 100.0
# Br8667_ant   99.4   0.6 100.0
# Br8667_mid   99.7   0.3 100.0
# Sum          91.6   8.4 100.0

#### QC plots ####
pdf(here("plots","03_build_sce","QC_outliers.pdf"), width = 21)
## Mito rate
plotColData(sce, x = "Sample", y="subsets_Mito_percent", colour_by="high_mito") +
  ggtitle("Mito Precent")+
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)

# ## low sum
plotColData(sce, x = "Sample", y="sum", colour_by="low_sum") +
  scale_y_log10() +
  ggtitle("Total count") +
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

# ## low detected
plotColData(sce, x = "Sample", y="detected", colour_by="low_detected") +
  scale_y_log10() +
  ggtitle("Detected features") +
  # geom_hline(yintercept = 500)+
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)

# Mito rate vs n detected features

plotColData(sce, x="detected", y="subsets_Mito_percent",
            colour_by="discard_auto", point_size=2.5, point_alpha=0.5) 

# Detected features vs total count

plotColData(sce, x="sum", y="detected",
            colour_by="discard_auto", point_size=2.5, point_alpha=0.5) 

dev.off()


#### Doublet detection ####
## To speed up, run on sample-level top-HVGs - just take top 1000 
set.seed(506)

colData(sce)$doubletScore <- NA

for(i in splitit(sce$Sample)){

  sce_temp <- sce[,i]
  ## To speed up, run on sample-level top-HVGs - just take top 1000 
  normd <- logNormCounts(sce_temp)
  geneVar <- modelGeneVar(normd)
  topHVGs <- getTopHVGs(geneVar, n=1000)
  
  dbl_dens <- computeDoubletDensity(normd, subset.row=topHVGs)
  colData(sce)$doubletScore[i] <- dbl_dens
  
}

summary(sce$doubletScore)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1439  0.4391  0.8215  1.0823 31.4607 

## Visualize doublet scores ##

dbl_df <- colData(sce) %>%
  as.data.frame() %>%
  select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
  ggplot(aes(x = reorder(Sample, doubletScore, FUN = median), y = doubletScore)) +
  geom_boxplot() +
  labs(x = "Sample") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  coord_flip() +
  my_theme

ggsave(dbl_box_plot, filename = here("plots","03_build_sce","doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
  ggplot(aes(x  = doubletScore)) +
  geom_density() +
  labs(x = "doublet score") +
  facet_grid(Sample ~ .) +
  theme_bw()

ggsave(dbl_density_plot, filename = here("plots","03_build_sce","doublet_scores_desnity.png"), height = 17)

dbl_df %>%
  group_by(Sample) %>%
  summarize(median = median(doubletScore),
            q95 = quantile(doubletScore, .95),
            drop = sum(doubletScore >= 5),
            drop_precent = 100*drop/n())

# Sample      median   q95  drop drop_precent
# <chr>        <dbl> <dbl> <int>        <dbl>
# 1 Br2720_mid   0.299  1.81    44       1.18  
# 2 Br2720_post  1.01   3.72    89       1.50  
# 3 Br2743_ant   0.400  1.79    23       0.719 
# 4 Br2743_mid   0.281  1.19    12       0.350 
# 5 Br3942_ant   0.376  1.42    51       0.814 
# 6 Br3942_mid   0.448  1.66     6       0.139 
# 7 Br6423_ant   0.376  1.98    30       0.750 
# 8 Br6423_post  0.301  1.82    23       0.566 
# 9 Br6432_ant   0.474  2.45    33       0.962 
# 10 Br6471_ant   0.355  1.62    35       0.966 
# 11 Br6471_mid   0.278  1.48    46       0.859 
# 12 Br6522_mid   1.62   4.71   145       3.59  
# 13 Br6522_post  1.85   5.72   352       8.20  
# 14 Br8325_ant   0.321  1.51    60       1.01  
# 15 Br8325_mid   0.250  1.19    34       0.681 
# 16 Br8492_mid   0.477  2.18    10       0.193 
# 17 Br8492_post  0.496  1.80     1       0.0335
# 18 Br8667_ant   0.453  2.23    33       0.568 
# 19 Br8667_mid   0.451  2.44    14       0.339 


table(sce$discard_auto, sce$doubletScore >= 5)
#       FALSE  TRUE
# FALSE 76581  1023
# TRUE   7134    18


## Save
save(sce, file = here::here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

## Save out sample info for easy access
pd <- colData(sce) %>% as.data.frame()

sample_info <- pd %>% 
  group_by(Sample, file_id, region, subject, round) %>%
  summarize(n = n(),
            n_high_mito = sum(high_mito),
            n_low_sum = sum(low_sum),
            n_low_detect = sum(low_detected),
            n_discard_auto = sum(discard_auto)
            )

write_csv(sample_info, file = here("processed-data", "03_build_sce","sample_info.csv"))

n_boxplot <- sample_info %>%
  ggplot(aes(x= round, y = n, color = round)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  ggrepel::geom_text_repel(aes(label = Sample), color = "black") +
  my_theme +
  theme(legend.position = "None")

ggsave(n_boxplot , filename = here("plots","03_build_sce", "droplet_qc" ,"n_nuclei_boxplot.png"))

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
