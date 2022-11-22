library(SpatialExperiment)
library(purrr)

# Load spe of all samples ------------------------------------------------------
path_nonIF_in <- paste(
    "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC",
    "processed-data", "rdata", "spe", "01_build_spe",
    "spe_filtered_final_with_clusters.Rdata", sep = "/"
)
load(path_nonIF_in)
# NOTE: data are load to the environment names ***spe***

## Number of samples
Sample_IDs <- spe$sample_id |> unique()


# Create sample-specific spe ---------------------------------------
if(!dir.exists("~/CCC_Visium/")) dir.create("~/CCC_Visium/", recursive = TRUE)
Sample_IDs |>
    purrr::map(.f = function(ID){
        folder_name <- paste0("~/CCC_Visium/", ID)
        # Create folder
        if(!dir.exists(folder_name))
            dir.create(folder_name, recursive = TRUE)
        file_path <- paste0(folder_name,"/spe.rds")
        # Create spe rds
        if(!file.exists(file_path)){
            spe[, spe$sample_id == ID] |>
                saveRDS(file = file_path)
        }
    })



