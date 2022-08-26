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

## plot dir
plot_dir <- here("plots", "05_explore_sce", "07_spatial_registration")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# load sc data
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## previously saved data
# load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity.RDS", verbose = TRUE)
load(here("processed-data", "05_explore_sce", "enrichment", "enrichment_cellType_hc.Rdata"))
# results_specificity
## Extract t-stats
specificity_stats <- results_specificity[, grep("^t_stat", colnames(results_specificity))]
colnames(specificity_stats) <- gsub("^t_stat_", "", colnames(specificity_stats))

dim(specificity_stats)
# [1] 5544   29

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

## plot
pdf(here(plot_dir, "spatial_registration_plot_hc_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_hc_top100, max = max(cor_hc_top100))
dev.off()

#### Excit Only ####
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/12_spatial_registration_sc/results_specificity_excit.RDS", verbose = TRUE)
## Fix head code bug for now
colnames(results_specificity_excit)[[11]] <- "t_stat_Excit_10"
specificity_stats_excit <- results_specificity_excit[, grep("^t_stat", colnames(results_specificity_excit))]
colnames(specificity_stats_excit) <- gsub("^t_stat_", "", colnames(specificity_stats_excit))

## HC anno vs. top 100
cor_excit_top100 <- layer_stat_cor(
    specificity_stats_excit,
    modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    top_n = 100
)

## Plot
pdf(here(plot_dir, "spatial_registration_plot_Excit_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_excit_top100, max = max(cor_excit_top100))
dev.off()

#### Azimuth registration ####
load(here("processed-data", "05_explore_sce", "enrichment", "enrichment_cellType_azimuth.Rdata"))

specificity_stats_azimuth <- results_specificity[, grep("^t_stat", colnames(results_specificity))]
colnames(specificity_stats_azimuth) <- gsub("^t_stat_", "", colnames(specificity_stats_azimuth))

## HC anno vs. top 100
cor_azimuth_top100 <- layer_stat_cor(
    specificity_stats_azimuth,
    modeling_results,
    model_type = names(modeling_results)[2],
    reverse = FALSE,
    top_n = 100
)

## Plot
pdf(here(plot_dir, "spatial_registration_plot_azimuth_v_manual_top100.pdf"))
layer_stat_cor_plot(cor_azimuth_top100, max = max(cor_azimuth_top100))
dev.off()

## Save correlations
save(cor_hc_top100, cor_excit_top100, cor_azimuth_top100, file = here("processed-data", "05_explore_sce", "07_spatial_registration", "layer_cor.Rdata"))
