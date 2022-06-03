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
load(here("processed-data", "03_build_sce","sce_MNN_round.Rdata"), verbose = TRUE)
load(here("processed-data", "03_build_sce", "clusters.Rdata"), verbose = TRUE)

# Assign as 'prelimCluster'
sce$prelimCluster <- factor(clusters)
plotReducedDim(sce, dimred="TSNE", colour_by="prelimCluster")

table(sce$prelimCluster)

## Many small clusters, some very large clusters
summary(as.numeric(table(sce$prelimCluster)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7.0    53.0   105.0   261.3   221.0 11945.0 

# Is sample driving this 'high-res' clustering at this level?
round_prelimClusters <-table(sce$prelimCluster, sce$round)
round_prelimClusters[which(rowSums(round_prelimClusters != 0) == 1),]
## 270 clusters from round1
length(round_prelimClusters[which(rowSums(round_prelimClusters != 0) == 1),])

sample_prelimClusters <- table(sce$prelimCluster, sce$Sample)
sample_prelimClusters[which(rowSums(sample_prelimClusters != 0) == 1),]
# Br6432_ant 2 exclusive, Br6471_ant 7 exclusive

# # rbind the ref.sampleInfo[.rev]
# ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

#### check doublet score for each prelim clust ####
clusIndexes = splitit(sce$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.02322  0.29750  0.52881  0.86885  0.86777 22.70270 

clust_medianDoublet_df <- stack(prelimCluster.medianDoublet) %>%
  rename(medianDoublet = values, prelimCluster = ind)

medianDoublet_histo <- ggplot(clust_medianDoublet_df, aes(x = medianDoublet)) +
  geom_histogram() +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed") +
  my_theme

ggsave(medianDoublet_histo, filename = here(plot_dir, "medianDoublet_histogram.png"))

## watch in clustering
(doublet_clust <- prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5])
# 30       167       233       257       272 
# 22.702696 14.723280  6.450884  5.065352 13.055584 

table(sce$prelimCluster)[names(doublet_clust)]
# 30 167 233 257 272 
# 39  19  17  32  22 


#### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization ####
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing

# Preliminary cluster index for pseudo-bulking
prelimCluster.PBcounts <- sapply(clusIndexes, function(ii){
  rowSums(assays(sce)$counts[ ,ii])
}
)

dim(prelimCluster.PBcounts)
# [1] 36601   297

corner(prelimCluster.PBcounts)

# And btw...
table(rowSums(prelimCluster.PBcounts)==0)
# FALSE  TRUE 
# 34987  1614 

# Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(prelimCluster.PBcounts)

# Normalize with these LSFs
geneExprs.temp <- t(apply(prelimCluster.PBcounts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

# Print for future reference
pdf(here(plot_dir, "dend.pdf"), height = 12)
    par(cex=0.6, font=2)
    plot(dend, main="test dend", horiz = TRUE)
    abline(v = 525, lty = 2)
dev.off()
    

# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

# Just for observation
par(cex=.6)
myplclust(tree.clusCollapsed, cex.main=2, cex.lab=1.5, cex=1.8)

dend %>%
  set("labels_cex", 0.8) %>%
  plot(horiz = TRUE)
abline(v = 525, lty = 2)

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=525)

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
pdf(here(".pdf", height = 9)
par(cex=0.6, font=2)
plot(dend, main="3x DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = 325, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                  merged=clust.treeCut[order.dendrogram(dend)])


# Assign as 'collapsedCluster'
sce$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce$prelimCluster, clusterRefTab.dlpfc$origClust)])
n_clusters <- length(levels(sce$collapsedCluster))
# Print some visualizations:
pdf("pdfs/revision/regionSpecific_DLPFC-n3_reducedDims-with-collapsedClusters_LAH2021.pdf")
plotReducedDim(sce, dimred="PCA_corrected", ncomponents=5, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce, colour_by="sampleID", point_alpha=0.5)
plotTSNE(sce, colour_by="protocol", point_alpha=0.5)
plotTSNE(sce, colour_by="collapsedCluster", point_alpha=0.5)
plotTSNE(sce, colour_by="sum", point_alpha=0.5)
plotTSNE(sce, colour_by="doubletScore", point_alpha=0.5)
# And some more informative UMAPs
plotUMAP(sce, colour_by="sampleID", point_alpha=0.5)
plotUMAP(sce, colour_by="collapsedCluster", point_alpha=0.5)
dev.off()

tail(table(sce$prelimCluster, sce$collapsedCluster),20)
