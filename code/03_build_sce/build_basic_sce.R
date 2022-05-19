# qrsh -l bluejay,mem_free-80G,h_vmem=80G

library("SingleCellExperiment")
library("SpatialExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")

## Define some info for the samples
tmp <- read.table(here("raw-data", "sample_libs_round0-5.tsv"), header = FALSE)

sample_info <- data.frame(
    file_id = tmp$V1,
    region_short = tolower(gsub("DLPFC_", "", tmp$V2)),
    subject = tmp$V3,
    round = tmp$V5
)
rm(tmp)
sample_info$region <- gsub("mid", "middle", sample_info$region_short)
sample_info$region[sample_info$region != "middle"] <- paste0(sample_info$region[sample_info$region != "middle"], "erior")
sample_info$Sample <- with(sample_info, paste0(subject, "_", region_short))
stopifnot(all(!duplicated(sample_info$Sample)))

sample_info$sample_path <- file.path(
    here::here("processed-data", "cellranger"),
    sample_info$file_id,
    "outs",
    "raw_feature_bc_matrix"
)
stopifnot(all(file.exists(sample_info$sample_path)))

## Get the rest of the donor info from the spatialDLPFC SPE object
## Load SPE data
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/code/deploy_app/spe_merged_final_nocounts.Rdata", verbose = TRUE)
m <- match(sample_info$subject, spe$subject)
sample_info$age <- spe$age[m]
sample_info$sex <- spe$sex[m]
sample_info$diagnosis <- spe$diagnosis[m]
rm(m, spe)

## Build basic SCE
Sys.time()
sce <- read10xCounts(
    sample_info$sample_path,
    sample_info$Sample,
    type = "sparse",
    col.names = TRUE
)
Sys.time()
# [1] "2021-11-30 10:25:04 EST"
# [1] "2021-11-30 10:40:51 EST"

## Use key similar to spe objects
sce$key <- paste0(sce$Barcode, "_", sce$Sample)

## Add the study design info
new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
## Fix order
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SPE object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce
# class: SingleCellExperiment
# dim: 36601 26605806
# metadata(1): Samples
# assays(1): counts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817 ENSG00000277196
# rowData names(6): source type ... gene_name gene_type
# colnames(26605806): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ... 19_TTTGTTGTCTTTGGCT-1
#   19_TTTGTTGTCTTTGTCG-1
# colData names(11): Sample Barcode ... sex diagnosis
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

## Note that no empty droplets have been filtered out yet!

## Save for later
dir.create(here::here("processed-data", "sce"), showWarnings = FALSE)
save(sce, file = here::here("processed-data", "sce", "sce_raw.Rdata"))

## Size in Gb
lobstr::obj_size(sce) / 1024^3
# 15.52167

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info  ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  hash: woman factory worker: medium-dark skin tone, woman and man holding hands, raising hands: dark skin tone
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
#  BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
#  BiocParallel           1.28.2   2021-11-25 [2] Bioconductor
#  Biostrings             2.62.0   2021-10-26 [2] Bioconductor
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
#  GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
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
#  lobstr               * 1.1.1    2019-07-02 [2] CRAN (R 4.1.0)
#  locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
#  magick                 2.7.3    2021-08-18 [2] CRAN (R 4.1.2)
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
#  restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
#  rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
#  Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
#  rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
#  rmote                  0.3.4    2021-11-02 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools              2.10.0   2021-10-26 [2] Bioconductor
#  rtracklayer            1.54.0   2021-10-26 [2] Bioconductor
#  S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  scuttle                1.4.0    2021-10-26 [1] Bioconductor
#  servr                  0.23     2021-08-11 [1] CRAN (R 4.1.2)
#  sessioninfo          * 1.2.1    2021-11-02 [2] CRAN (R 4.1.2)
#  SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
#  sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
#  SpatialExperiment    * 1.4.0    2021-10-26 [1] Bioconductor
#  SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
#  tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  xfun                   0.28     2021-11-04 [2] CRAN (R 4.1.2)
#  XML                    3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
#  XVector                0.34.0   2021-10-26 [2] Bioconductor
#  yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.1.x
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
