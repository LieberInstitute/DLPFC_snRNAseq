library("SingleCellExperiment")
library("SpatialExperiment")
library("DropletUtils")
library("here")
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
sample_info$sample_id <- with(sample_info, paste0(subject, "_", region_short))
stopifnot(all(!duplicated(sample_info$sample_id)))

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
  sample_info$sample_id,
  type = "sparse",
  col.names = TRUE
)
Sys.time()
# [1] "2021-11-29 16:41:53 EST"
# 

## Add the study design info
key_original <- colnames(sce)
new_col <- merge(colData(sce), sample_info)
## Fix order
new_col <- new_col[match(key_original, rownames(new_col)), ]
stopifnot(identical(rownames(new_col), key_original))
rownames(new_col) <- rownames(colData(spe))
colData(sce) <- new_col[, -which(colnames(new_col) == "sample_path")]

## Save for later
dir.create(here::here("processed-data", "sce"), showWarnings = FALSE)
save(sce, file = here::here("processed-data", "sce", "sce_raw.Rdata"))

## Size in Gb
lobstr::obj_size(sce) / 1024^3

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

