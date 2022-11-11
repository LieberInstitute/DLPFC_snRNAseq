library(tidyverse)
library(SingleCellExperiment)


# Load Batch Effect Correct Data
load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/03_build_sce/sce_harmony_Sample.Rdata")
# sce <- readRDS(file = "processed-data/03_build_sce/DLPFC_sce.rds")
col_data <- colData(sce) |> as.data.frame()


# Sample Summary Statistics
col_data |> group_by(subject) |> summarize(n = n_distinct(region))


# Sample Br2720
# Why when I do this subsetting, we have loading all these packages again?
# Which object is doing this annoying thing!
sce_Br2720 <- sce[, sce$subject=="Br2720"]


library(RSoptSC)

