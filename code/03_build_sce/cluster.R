library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran") 
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")

## Best normalization result
load(here("processed-data", "03_build_sce","sce_MNN_subject.Rdata"))

message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k=10, use.dimred="PCA")

message("running walktrap - ", Sys.time())
clusters <- igraph::cluster_walktrap(snn.gr)$membership

table(clusters)
message("saving data - ", Sys.time())
save(clusters, here("processed-data", "03_build_sce", "clusters.Rdata"), verbose = TRUE)

## Save final sce w/ annotations
# save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('cluster', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript cluster.R subject")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
