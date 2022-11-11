library(tidyverse)
library(SingleCellExperiment)


# Load Batch Effect Correct Data
# load(file = "processed-data/03_build_sce/sce_harmony_Sample.Rdata")
sce <- readRDS(file = "processed-data/03_build_sce/DLPFC_sce.rds")
col_data <- colData(sce) |> as.data.frame()


# Sample Summary Statistics
col_data |> group_by(subject) |> summarize(n = n_distinct(region))


# Sample Br2720
sce_Br2720 <- sce[, sce$subject=="Br2720"]


