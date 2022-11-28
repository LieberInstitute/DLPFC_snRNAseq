library(tidyverse)
library(SingleCellExperiment)


# Receive Arguments from Command Line
args <- (commandArgs(TRUE))

if (length(args) == 0) {
    print("No arguments supplied.")
} else {
    for (i in 1:length(args)) {
        eval(parse(text = args[[i]]))
    }
}

## NOTE: reserved for when command line arguments fails
# crn_Positions <- c("Anterior", "Middle", "Posterior")
# crn_sec <- crn_Positions[1]

# Create section specific folders to contain results
CCC_res_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/07_ccc/"
if (!dir.exists(CCC_res_path)) dir.create(CCC_res_path, recursive = TRUE)
fdl_path <- paste0(CCC_res_path, crn_sec, "/")
if (!dir.exists(fdl_path)) dir.create(fdl_path, recursive = TRUE)

# Data Prep ---------------------------------------------------------------

# Load Proprocessed Data
## NOTE: Set working directory is necessary in order to use the rds data because of the .h5 file
## https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1668537164053579?thread_ts=1668530634.509189&cid=C01EA7VDJNT
setwd("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/")
sce <- readRDS("se.rds")

# Summary stat for coronary sections
colData(sce) |>
    as.data.frame() |>
    group_by(Sample) |>
    summarize(n = n())


## Check celltype columns definition
# cData <- colData(sce) |> data.frame()
#
# colData(sce) |> data.frame() |>
#     select(starts_with("cellType")) |>
#     summarize(across(.fns = n_distinct))
# xtabs(~ cellType_broad_k + cellType_layer, cData, addNA = TRUE)

# Subset cells from each coronal section
# R2 Code
# sce_crn <- sce[, sce$Position==crn_sec]
# R3 Code
sce_crn <- sce[, sce$Sample == crn_sec]
# if(ncol(sce_crn) <=0) stop

# QC: Remove bad cells----------------------------------------------------------------------
# Following Louise's suggstion, remove drop cells
# https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1668530673111519?thread_ts=1668530634.509189&cid=C01EA7VDJNT
sce_crn <- sce_crn[, sce_crn$cellType_hc != "drop"]
sce_crn$cellType_hc <- droplevels(sce_crn$cellType_hc)




## See if batch effect is large
# library(scater)
# ggsave(paste0( fdl_path,"tsne.png"),
#        plot = plotReducedDim(sce_crn, dimred="TSNE", colour_by="BrNum"))

## Empty Cells/Genes
# sum(rowSums(logcounts(sce_crn))==0)
# sum(colSums(logcounts(sce_crn))==0)


# HVG
# TODO: Decides if we want to prioritize some genes to reduce computation load
# TODO: If so, change the sce object names
# library(scran)
# dec_sce <- modelGeneVar(sce_Br2720)
# hvg_names <- getTopHVGs(dec_sce, n=100)

# LIANA Analysis ----------------------------------------------------------
library(liana)
# Prep
## Assigning cellType to sce label (Required by LIANA)
## Remove NA cellType
sce_crn <- sce_crn[, !is.na(colData(sce_crn)$cellType_layer)]
colLabels(sce_crn) <- colData(sce_crn)$cellType_layer

# TEST CODE: small scale
# sce_crn <- sce_crn[, 1:1000]


# Run liana
liana_test <- liana_wrap(sce_crn)
saveRDS(liana_test,
    file = paste0(fdl_path, "liana_test.rds")
)


# If segfault happens, remove that list in the list, e.g.
# liana_test <- liana_test
# liana_test$cellphonedb <- NULL

# Consensus Ranks
liana_res <- liana_test %>%
    liana_aggregate()
saveRDS(liana_res,
    file = paste0(fdl_path, "liana_consensus.rds")
)

sessionInfo()
