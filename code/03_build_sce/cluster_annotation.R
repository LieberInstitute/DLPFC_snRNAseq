library("SingleCellExperiment")
# library(EnsDb.Hsapiens.v86)
library("scater")
# library(scran)
# library(batchelor)
# library(DropletUtils)
library("jaffelab")
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")

## Load sce object
load(here("processed-data", "03_build_sce","sce_MNN_subject.Rdata"), verbose = TRUE)
load(here("processed-data", "03_build_sce", "clusters.Rdata"), verbose = TRUE)


# Assign as 'prelimCluster'
sce.dlpfc$prelimCluster <- factor(clusters.k20)
plotReducedDim(sce.dlpfc, dimred="TSNE", colour_by="prelimCluster")

# Is sample driving this 'high-res' clustering at this level?
(sample_prelimClusters <- table(sce.dlpfc$prelimCluster, sce.dlpfc$sampleID))  # (a little bit, but is typical)
sample_prelimClusters[which(rowSums(sample_prelimClusters == 0) == 2),]
# 39 - only 4 samples all from Br5207

# rbind the ref.sampleInfo[.rev]
ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

## check doublet score for each prelim clust
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce.dlpfc$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.01059  0.07083  0.14823  0.53264  0.30064 14.79144

hist(prelimCluster.medianDoublet)

## watch in clustering
prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5]
# 19       32       73
# 14.79144  7.98099 10.20462

table(sce.dlpfc$prelimCluster)[c(19, 32, 73)]
# 19 32 73
# 27 32  8

# Save for now
save(sce.dlpfc, chosen.hvgs.dlpfc, pc.choice.dlpfc, ref.sampleInfo,
     file="/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda")


### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing

# Preliminary cluster index for pseudo-bulking
clusIndexes = splitit(sce.dlpfc$prelimCluster)
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce.dlpfc)$counts[ ,ii])
}
)

# And btw...
table(rowSums(prelimCluster.PBcounts)==0)
# FALSE  TRUE
# 29310  4228

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>%
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 325, lty = 2)

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=325)

table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
## Cutting at 250 looks good for the main neuronal branch, but a lot of glial
#    prelim clusters are dropped off (0's)

# # Cut at 400 for broad glia branch (will manually merge remaining dropped off)
# glia.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
#                               minClusterSize=2, deepSplit=1, cutHeight=400)
# unname(glia.treeCut[order.dendrogram(dend)])

# Take those and re-assign to the first assignments

# clust <- clust.treeCut[order.dendrogram(dend)]
# clust2 <- name_zeros(clust, list(c(1,2), c(106,107)))
# unname(clust2)

# Add new labels to those prelimClusters cut off
clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max(clust.treeCut)+c(1, 6, 2, 3, 4, 5, 5)

# 'Re-write', since there are missing numbers
# clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust2))
clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))

## Define color pallet
cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])
labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf("pdfs/revision/regionSpecific_DLPFC-n3_HC-prelimCluster-relationships_LAH2021.pdf", height = 9)
par(cex=0.6, font=2)
plot(dend, main="3x DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                  merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce.dlpfc$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce.dlpfc$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce.dlpfc$collapsedCluster))
# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
plotReducedDim(sce.dlpfc, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="sum", point_alpha=0.5)
plotTSNE(sce.dlpfc, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce.dlpfc, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce.dlpfc, colour_by="collapsedCluster", point_alpha=0.5)
dev.off()

tail(table(sce.dlpfc$prelimCluster, sce.dlpfc$collapsedCluster),20)
