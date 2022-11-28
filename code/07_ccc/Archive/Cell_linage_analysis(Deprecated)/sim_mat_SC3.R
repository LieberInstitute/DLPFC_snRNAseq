library(tidyverse)
library(SingleCellExperiment)


# Data Prep ---------------------------------------------------------------

# Load Batch Effect Correct Data
# load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/03_build_sce/sce_harmony_Sample.Rdata")
load(file = "processed-data/03_build_sce/sce_harmony_Sample.Rdata")
# sce <- readRDS(file = "processed-data/03_build_sce/DLPFC_sce.rds")
# col_data <- colData(sce) |> as.data.frame()
#
#
# # Sample Summary Statistics
# col_data |> group_by(subject) |> summarize(n = n_distinct(region))


# Sample Br2720
# Why when I do this subsetting, we have loading all these packages again?
# Which object is doing this annoying thing!
sce_Br2720 <- sce[, sce$subject == "Br2720"]



## Library Size Normalization
library(scater)
set.seed(1000)
lib_Br2720 <- librarySizeFactors(sce_Br2720)
sizeFactors(sce_Br2720) <- lib_Br2720
sce_Br2720 <- logNormCounts(sce_Br2720)

## HVG
library(scran)
dec_sce <- modelGeneVar(sce_Br2720)
# Examine Mean-variance relationship
# fit.pbmc <- metadata(dec_sce)
# pdf("~/mean_var.pdf")
# plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
#      ylab="Variance of log-expression")
# curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)
# dev.off()
## The mean variance relationship looks normal that is abnormal

# Highly Variable Genes
# TODO: edit this to a much larger number
hvg_names <- getTopHVGs(dec_sce, n = 100)

# TODO: remove the filter on number of cells
hvg_sce_Br2720 <- sce_Br2720[hvg_names, 1:500]


# Calculate the similarity matrix with SC3 ----


## Prepare default rowData
rowData(hvg_sce_Br2720)$feature_symbol <- rowData(hvg_sce_Br2720)$gene_name

convert_sparse_matrix <- function(sce) {
    # TODO: change hard coding conversion to else
    counts(sce) <- as.matrix(counts(sce))
    logcounts(sce) <- as.matrix(logcounts(sce))
    return(sce)
}

# Remove
sum((logcounts(hvg_sce_Br2720) |> as.matrix() |> colSums()) == 0)
sum((logcounts(hvg_sce_Br2720) |> as.matrix() |> rowSums()) == 0)

hvg_sce_Br2720 <- hvg_sce_Br2720[, (logcounts(hvg_sce_Br2720) |> as.matrix() |> colSums()) != 0]

library(SC3)
# sc3_sce <- sc3(convert_sparse_matrix(hvg_sce_Br2720), ks = 3, biology = TRUE, gene_filter = TRUE, n_core = 4)
set.seed(1)
sc3_sce <- sc3_prepare(convert_sparse_matrix(hvg_sce_Br2720), gene_filter = TRUE, n_core = 4)

# TODO:
# sc3_estimate_k

sc3_sce <- sc3_calc_dists(sc3_sce)
sc3_sce <- sc3_calc_transfs(sc3_sce)
