library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")
library("dplyr")

# Plotting set up 
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","cluster")

if(!dir.exists(plot_dir)) dir.create(plot_dir)

## Load sce object
load(here("processed-data", "03_build_sce","sce_harmony_Sample.Rdata"), verbose = TRUE)
load(here("processed-data", "03_build_sce", "clusters.Rdata"), verbose = TRUE)

# Assign as 'prelimCluster'
sce$prelimCluster <- factor(clusters)
length(levels(sce$prelimCluster))
# [1] 296
table(sce$prelimCluster)

clusIndexes = splitit(sce$prelimCluster)

message("Pseudobulk - ", Sys.time())
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce)$counts[ ,ii])
}
)

dim(prelimCluster.PBcounts)
# [1] 36601   297

# And btw...
table(rowSums(prelimCluster.PBcounts)==0)

message("Get Lib Size Factors - ", Sys.time())
# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
message("Normalize - ", Sys.time())
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
message("Cluster Again - ", Sys.time())
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Print for future reference
pdf(here(plot_dir, "dend.pdf"), height = 12)
    par(cex=0.6, font=2)
    plot(dend, main="hierarchical cluster dend", horiz = TRUE)
    # abline(v = 525, lty = 2)
dev.off()

## Save data
save(dend, tree.clusCollapsed, dist.clusCollapsed, file = here("processed-data", "03_build_sce", "HC_dend.Rdata"))

# sgejobs::job_single('hierarchical_cluster', create_shell = TRUE, memory = '100G', command = "Rscript hierarchical_cluster.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
