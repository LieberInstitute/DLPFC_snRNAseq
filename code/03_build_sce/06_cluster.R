library("SingleCellExperiment")
library("jaffelab")
library("scater")
library("scran") 
library("here")
library("sessioninfo")

## Best normalization result
load(here("processed-data", "03_build_sce","sce_harmony_Sample.Rdata"), verbose = TRUE)

message("running buildSNNGraph - ", Sys.time())
snn.gr <- buildSNNGraph(sce, k=20, use.dimred="HARMONY")

message("running walktrap - ", Sys.time())
clusters <- igraph::cluster_walktrap(snn.gr)$membership

table(clusters)
message("saving data - ", Sys.time())
save(clusters, file = here("processed-data", "03_build_sce", "clusters.Rdata"))

## Save final sce w/ annotations
# save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('06_cluster', create_shell = TRUE, queue= 'bluejay', memory = '100G', command = "Rscript 06_cluster.R subject")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
