library(tidyverse)
library(readr)
library(SpatialExperiment)
# library(unglue)

deconvo_res_path <- paste0(
    "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/",
    "processed-data/spot_deconvo/05-shared_utilities/nonIF/",
    "results_raw_layer.csv"
)


tangram_res <- read.table(deconvo_res_path,sep = ",", header  = TRUE) |>
    filter(deconvo_tool == "tangram") #|>
    #mutate(key = paste(barcode, sample_id, sep = "_"))

# Load spe of all samples ------------------------------------------------------

# The clustering raw results have inconsistent naming convention,
# Cause problem when merging
# https://github.com/LieberInstitute/spatialDLPFC/issues/133
# cluster_fld_path <- paste(
#     "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC",
#     "processed-data", "rdata", "spe", "clustering_results", sep = "/"
# )
#
# sp9_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp9 <- read.csv(sp9_path, header = TRUE)
#
#
# sp9_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp9 <- read.csv(sp9_path, header = TRUE) |>
#     dplyr::rename(sp9 = cluster)
#
# # Code to debug sp9 cluster loading
# # tmp <- sp9 |> unglue_unnest(col = key,
# #                      patterns = "{spot_id}_{Br_num}_{pos}")
# #
# # tmp |> filter(Br_num == "Br2720")
#
# sp16_path <- paste(cluster_fld_path,
#                   "bayesSpace_harmony_9",
#                   "clusters.csv", sep = "/")
# sp16 <- read.csv(sp16_path, header = TRUE)|>
#     dplyr::rename(sp16 = cluster)

load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/01_build_spe/spe_filtered_final_with_clusters.Rdata")

# tmp <- colData(spe) |> data.frame() |>
#     dplyr::select(key:array_col,
#                   sum_umi:col,
#                   bayesSpace_harmony_9, bayesSpace_harmony_16)

# rm(spe)

fnl_dat <- tangram_res |>
    full_join(
        # tmp,
    colData(spe) |> data.frame() |>
    dplyr::select(key:array_col,
                  sum_umi:col,
                  bayesSpace_harmony_9, bayesSpace_harmony_16),
    by = "sample_id"
    )

saveRDS(fnl_dat,
        file = "~/cell_comp_full_dat.rds")




