
library("SingleCellExperiment")
library("scater")
library("here")
library("sessioninfo")

load(here("processed-data","sce","sce-DLPFC.Rdata"), verbose = TRUE)

#### Find Markers ####
(small_clusters <- cluster_tab[cluster_tab < 20])
# mbk05 mbk18 mbk24 mbk29 
# 11     3     2     7

sce <- sce[,!(sce$kmeans %in% names(small_clusters))]
dim(sce)
# [1] 36601 77581

colLabels(sce) <- sce$kmeans
markers <- scran::findMarkers(sce, pval.type="all", direction="up")

save(markers, file = here("processed-data", "03_build_sce","kmeans_29_markers.Rdata"))
