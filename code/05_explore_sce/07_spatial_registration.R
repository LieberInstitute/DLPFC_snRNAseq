library("SpatialExperiment")
library("spatialLIBD")
library("HDF5Array")
# library("rafalib")
# library("scuttle")
# library(limma)
# library("RColorBrewer")
# library(lattice)
# library(edgeR)
library("here")
library("sessioninfo")

source("utils.R")

# load sc data
# load(file = "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC.Rdata")
sce <- loadHDF5SummarizedExperiment(here("processed-data","sce","sce_DLPFC_annotated"))

# var_oi <- "cellType_hc"
# covars <- c("Position", "age", "sex")

rownames(sce) <- rowData(sce)$gene_id # have to make row names of object the ensembl id instead of gene names
colData(sce)$Position <- as.factor(colData(sce)$Position)
colData(sce)$age <- as.numeric(colData(sce)$age)
colData(sce)$sex <- as.factor(colData(sce)$sex)
colnames(colData(sce))[1] <- "sample_id"

colData(sce)

spe_pseudo <- scuttle::aggregateAcrossCells(
  sce,
  DataFrame(
    BayesSpace = colData(sce)[["cellType_hc"]],
    sample_id = sce$sample_id
  )
)

results_specificity <- computeEnrichment(sce, var_oi = "cellType_hc", covars = c("Position", "age", "sex"))

# save(results_specificity, file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity.RDS")

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity.RDS")

specificity_stats <- results_specificity[, grep("^t_stat", colnames(results_specificity))]
colnames(specificity_stats) <- gsub("^t_stat_", "", colnames(specificity_stats))

dim(specificity_stats)
# [1] 5544   28

colnames(specificity_stats) # missing Excit_07
# [1] "Astro"         "Endo.Mural_01" "Endo.Mural_02" "Excit_01"      "Excit_02"      "Excit_03"      "Excit_04"     
# [8] "Excit_05"      "Excit_06"      "Excit_08"      "Excit_09"      "Excit_10"      "Excit_11"      "Excit_12"     
# [15] "Excit_13"      "Excit_14"      "Excit_15"      "Inhib_01"      "Inhib_02"      "Inhib_03"      "Inhib_04"     
# [22] "Inhib_05"      "Inhib_06"      "Micro"         "Oligo_01"      "Oligo_02"      "Oligo_03"      "OPC"        

# vs manual annotations
modeling_results <- fetch_data(type = "modeling_results")
cor <- layer_stat_cor(
  specificity_stats,
  modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE,
  top_n = NULL
)

cor
save(cor, file = here("processed-data","05_explore_sce","07_sptatial_registration","layer_cor.Rdata"))

pdf("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/plots/12_spatial_registration_sc/spatial_registration_plot_sc_v_manual.pdf")
layer_stat_cor_plot(cor, max = 0.7)
dev.off()

## Azimuth registration ##
results_specificity <- computeEnrichment(sce, var_oi = "cellType_azimuth", covars = c("Position", "age", "sex"))

