library("SingleCellExperiment")
library("scater")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

#### Find Markers ####
## Drop NA in layer anno
sce <- sce[, !is.na(sce$cellType_layer)]

table(sce$cellType_layer)
# Astro Endo.Mural_01 Endo.Mural_02      Excit_01      Excit_02      Excit_03      Excit_04      Excit_05
# 3979           446          1711          7927          2487          1309          2171          2532
# Excit_06      Excit_07      Excit_08      Excit_09      Excit_10      Excit_11      Excit_12      Excit_13
# 329           334          1463          2561          1079           482           420          1567
# Excit_14      Excit_15      Inhib_01      Inhib_02      Inhib_03      Inhib_04      Inhib_05      Inhib_06
# 82            66          5366          1267          1310           565          1192          1367
# Micro      Oligo_01      Oligo_02      Oligo_03           OPC
# 1601         23025          4732          4294          1940

colLabels(sce) <- sce$cellType_layer

message("Running 1vALL findMarkers - ", Sys.time())
markers_1vALL <- scran::findMarkers(sce, pval.type = "all", direction = "up")


message("Running Tran 1vAll - ", Sys.time())
markers_1vALL_enrich <- findMarkers_1vAll(sce, assay_name = "logcounts", cellType_col = "cellType_layer", mod = "~Sample")

## Run mean Ratio markers
message("Running Mean Ratio - ", Sys.time())
markers_mean_ratio <- get_mean_ratio2(sce, assay_name = "logcounts", cellType_col = "cellType_layer")


save(markers_1vALL, markers_1vALL_enrich, markers_mean_ratio, file = here("processed-data", "03_build_sce", "cell_type_markers_layer.Rdata"))


# sgejobs::job_single('09_find_markers_layer', create_shell = TRUE, memory = '20G', command = "Rscript 09_find_markers_layer.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
