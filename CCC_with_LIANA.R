library(tidyverse)
library(SingleCellExperiment)


# Data Prep ---------------------------------------------------------------

# Load Batch Effect Correct Data
# load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/03_build_sce/sce_harmony_Sample.Rdata")
# load(file = "processed-data/03_build_sce/sce_harmony_Sample.Rdata")

sce <- readRDS("se.rds")

# col_data <- colData(sce) |> as.data.frame()
#
#
# # Sample Summary Statistics
# col_data |> group_by(BrNum) |> summarize(n = n_distinct(pos))



# QC ----------------------------------------------------------------------
# TODO: ask Louise and remove bad cells
# https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1668530673111519?thread_ts=1668530634.509189&cid=C01EA7VDJNT



# Sample Br2720
# TODO: change this to coronal sections
sce_Br2720 <- sce[, sce$BrNum=="Br2720"]



## Library Size Normalization
# library(scater)
# set.seed(1000)
# lib_Br2720 <- librarySizeFactors(sce_Br2720)
# sizeFactors(sce_Br2720)<-lib_Br2720
# sce_Br2720 <- logNormCounts(sce_Br2720)

## HVG
# library(scran)
# dec_sce <- modelGeneVar(sce_Br2720)
# # Examine Mean-variance relationship
# # fit.pbmc <- metadata(dec_sce)
# # pdf("~/mean_var.pdf")
# # plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
# #      ylab="Variance of log-expression")
# # curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)
# # dev.off()
# ## The mean variance relationship looks normal that is abnormal
#
# # Highly Variable Genes
# #TODO: edit this to a much larger number
# hvg_names <- getTopHVGs(dec_sce, n=100)
#
# #TODO: remove the filter on number of cells
# hvg_sce_Br2720 <- sce_Br2720[hvg_names,1:500]

# LIANA Analysis ----------------------------------------------------------

library(liana)
ccc_sce <- sce_Br2720[,1:500]

## Assigning cellType to sce label (Required by LIANA)
## TODO: (1) Move this label assigning step to the very begining of the process
## TODO: (2) decide what the cell type this analysis should be doing
colLabels(ccc_sce) <- colData(ccc_sce)$cellType_broad_hc

# Run liana
liana_test <- liana_wrap(ccc_sce)

# If segfault happens, remove that list in the list.
 tmp <- liana_test
tmp$cellphonedb <- NULL

# Consensus Ranks
tmp <- tmp %>%
    liana_aggregate()


library(circlize)
pdf(file = "~/CCC_toy.pdf")
chord_freq(tmp#,
                # source_groups = c("CD8 T", "NK", "B"),
                # target_groups = c("CD8 T", "NK", "B")
                )
dev.off()
