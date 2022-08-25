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
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

## previously saved data
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity.RDS", verbose = TRUE)
# results_specificity
head(results_specificity[,grep("Excit_07",colnames(results_specificity))])
#                 p_value_Excit_07 fdr_Excit_07 # No T-stats for Excit_07 ???
# ENSG00000000419        0.8128746    0.9210253
# ENSG00000001084        0.1046406    0.3697436
# ENSG00000001461        0.1951751    0.4999134
# ENSG00000001629        0.4219413    0.7080032
# ENSG00000001630        0.9801410    0.9959250
# ENSG00000001631        0.4906153    0.7499232

## Extract t-stats
specificity_stats <- results_specificity[, grep("^t_stat", colnames(results_specificity))]
colnames(specificity_stats) <- gsub("^t_stat_", "", colnames(specificity_stats))

dim(specificity_stats)
# [1] 5544   28

colnames(specificity_stats) # missing Excit_07
# [1] "Astro"         "Endo.Mural_01" "Endo.Mural_02" "Excit_01"      "Excit_02"      "Excit_03"      "Excit_04"     
# [8] "Excit_05"      "Excit_06"      "Excit_08"      "Excit_09"      "Excit_10"      "Excit_11"      "Excit_12"     
# [15] "Excit_13"      "Excit_14"      "Excit_15"      "Inhib_01"      "Inhib_02"      "Inhib_03"      "Inhib_04"     
# [22] "Inhib_05"      "Inhib_06"      "Micro"         "Oligo_01"      "Oligo_02"      "Oligo_03"      "OPC"        

#### cell types vs manual annotations ####
modeling_results <- fetch_data(type = "modeling_results")

## HC anno vs. top 100
cor_hc_top100 <- layer_stat_cor(
  specificity_stats,
  modeling_results,
  model_type = names(modeling_results)[2],
  reverse = FALSE,
  top_n = 100
)

save(cor_hc_top100 , file = here("processed-data","05_explore_sce","07_sptatial_registration","layer_cor.Rdata"))

## plot dir
plot_dir <- here("plots", "05_explore_sce",  "07_spatial_registration")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

## plot 
pdf(here(plot_dir,"spatial_registration_plot_hc_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_hc_top100, max = max(cor_hc_top100))
dev.off()

## Azimuth registration ##
results_specificity <- computeEnrichment(sce, var_oi = "cellType_azimuth", covars = c("Position", "age", "sex"))

