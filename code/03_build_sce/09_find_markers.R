library("SingleCellExperiment")
library("scater")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

#### Find Markers ####
table(sce$cellType_hc)
# Astro Endo.Mural_01 Endo.Mural_02      Excit_01      Excit_02      Excit_03      Excit_04      Excit_05
# 3979           446          1711          7927          2487          1309          2171          2532
# Excit_06      Excit_07      Excit_08      Excit_09      Excit_10      Excit_11      Excit_12      Excit_13
# 329           334          1463          2561          1079           482           420          1567
# Excit_14      Excit_15      Inhib_01      Inhib_02      Inhib_03      Inhib_04      Inhib_05      Inhib_06
# 82            66          5366          1267          1310           565          1192          1367
# Micro      Oligo_01      Oligo_02      Oligo_03           OPC
# 1601         23025          4732          4294          1940

colLabels(sce) <- sce$cellType_hc

message("Running 1vALL findMarkers")
Sys.time()
markers_1vALL <- scran::findMarkers(sce, pval.type = "all", direction = "up")
message("Done - ", Sys.time())
Sys.time()

## Run mean Ratio markers

markers_mean_ratio <- get_mean_ratio2(sce, cellType_col = "cellType_hc")
markers_mean_ratio_broad <- get_mean_ratio2(sce, cellType_col = "cellType_broad_hc")

save(markers_1vALL, markers_mean_ratio, markers_mean_ratio_broad, file = here("processed-data", "03_build_sce", "cell_type_markers.Rdata"))


# sgejobs::job_single('09_find_markers', create_shell = TRUE, queue= 'bluejay', memory = '20G', command = "Rscript 09_find_markers.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
