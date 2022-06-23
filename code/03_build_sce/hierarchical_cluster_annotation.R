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

## Many small clusters, some very large clusters
summary(as.numeric(table(sce$prelimCluster)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.0    54.0   109.0   262.2   222.2  8757.0

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
(doublet_clust <- prelimCluster.medianDoublet[prelimCluster.medianDoublet > 2])
# 25      156      166      199      247      258      264 
# 3.024884 8.395549 2.338546 2.330919 2.034900 2.159430 2.094400 
(doublet_clust <- prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5])
# 156 
# 8.395549 

summary(sce$doubletScore[sce$prelimCluster == 156])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3844  1.6301  8.3955  7.9715 13.7577 18.0338 

table(sce$prelimCluster)[names(doublet_clust)]
# 156 
# 68 


#### Step 2: Hierarchical clustering of pseudo-bulked ("PB'd") counts with most robust normalization ####
#         (as determined in: 'side-Rscript_testingStep2_HC-normalizn-approaches_wAmygData_MNTJan2020.R')
#           ** That is, to pseudo-bulk (aka 'cluster-bulk') on raw counts, on all [non-zero] genes,
#              normalize with `librarySizeFactors()`, log2-transform, then perform HC'ing


prelimCluster.PBcounts <- aggregateAcrossCells(sce,  DataFrame(
  prelimCluster = sce$prelimCluster
))

table(rowSums(assays(prelimCluster.PBcounts)$counts) == 0)
# # FALSE  TRUE
# # 34987  1614
#
# # Compute LSFs at this level
sizeFactors.PB.all  <- librarySizeFactors(assays(prelimCluster.PBcounts)$counts)
#
# # Normalize with these LSFs
geneExprs.temp <- t(apply(assays(prelimCluster.PBcounts)$counts, 1, function(x) {log2(x/sizeFactors.PB.all + 1)}))

## Perform hierarchical clustering
dist.clusCollapsed <- dist(t(geneExprs.temp))
tree.clusCollapsed <- hclust(dist.clusCollapsed, "ward.D2")

dend <- as.dendrogram(tree.clusCollapsed, hang=0.2)

## save HC data
save(dend, tree.clusCollapsed, dist.clusCollapsed, file = here("processed-data", "03_build_sce", "HC_dend.Rdata"))
# load(here("processed-data", "03_build_sce", "HC_dend.Rdata"), verbose = TRUE)

## pick k clusters?
# clus = cutree(dend, 19)
# cluster_colors <- c(c0 = "black", c(DeconvoBuddies::create_cell_colors(cell_types = paste0("c",1:19), pallet = "gg")))
# clus[156] <- 0
# 
# cluster_colors[clus + 1]
# 
# library(ape)
# 
# # Print for future reference
# pdf(here(plot_dir, "dend_aac.pdf"), height = 14)
#     par(cex=0.3)
#     plot(dend,main = "DLPFC prelim clusters", horiz = TRUE)
#     abline(v = 525, lty = 2)
# dev.off()


# pdf(here(plot_dir, "dend_test.pdf"), height = 14)
# par(cex=0.5)
# plot(as.phylo(dend),main = "DLPFC prelim clusters", tip.color = cluster_colors[clus +1])
# abline(v = 500, lty = 2)
# dev.off()


# labels(dend)[grep(c("19|32|73"),labels(dend))] <- paste0(labels(dend)[grep(c("19|32|73"),labels(dend))], "*")

## Try different cut heights

h_list <- seq(100, 2000, 25)

h_n_clust <- sapply(h_list, function(h) {
  clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                                 minClusterSize=2, deepSplit=1, cutHeight=h)
  return(length(table(clust.treeCut)))
})

ch_df <- data.frame(cut_height = h_list, n_clust = h_n_clust)

cut_height_scatter <- ggplot(ch_df, aes(x = cut_height, y = n_clust)) +
  geom_point() +
  my_theme +
  scale_x_reverse() +
  geom_vline(xintercept = 525, linetype = "dashed", color = "red")

ggsave(cut_height_scatter, filename = here(plot_dir, "cut_height_scatter.png"))

## Cut tree at chosen height

chosen_cut_height <- 800

clust.treeCut <- cutreeDynamic(tree.clusCollapsed, distM=as.matrix(dist.clusCollapsed),
                               minClusterSize=2, deepSplit=1, cutHeight=chosen_cut_height)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
unname(clust.treeCut[order.dendrogram(dend)]) == 0


# Add new labels to those prelimClusters cut off
## just define as a max cluster for now
if(any(clust.treeCut[order.dendrogram(dend)] == 0)){
  max_clust <- max(clust.treeCut) 
  clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)]==0)] <- max_clust + 1
  
  # 'Re-write', since there are missing numbers
  clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))
  cluster_colors[as.character(max_clust + 1)] <- "black"
  
}

## Define color pallet
# cluster_colors <- unique(tableau20[clust.treeCut[order.dendrogram(dend)]])
# names(cluster_colors) <- unique(clust.treeCut[order.dendrogram(dend)])

## too many groups for tableau (>20)
names(clust.treeCut) <- paste0("HC", numform::f_pad_zero(names(clust.treeCut)))
cluster_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(names(clust.treeCut))), pallet = "gg")


labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf(here(plot_dir, "dend_cut525.pdf"), height = 9)
par(cex=0.6, font=2)
plot(dend, main="DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = chosen_cut_height, lty = 2)
dev.off()


# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(origClust=order.dendrogram(dend),
                                  merged=names(clust.treeCut)[order.dendrogram(dend)])
head(clusterRefTab.dlpfc)

# Assign as 'collapsedCluster'
sce$collapsedCluster <- factor(clusterRefTab.dlpfc$merged[match(sce$prelimCluster, clusterRefTab.dlpfc$origClust)])

table(sce$collapsedCluster)

n_clusters <- length(levels(sce$collapsedCluster))
# Print some visualizations:
tail(table(sce$prelimCluster, sce$collapsedCluster),20)


## Plot clusters in TSNE
TSNE_clusters <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=collapsedCluster)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cluster_colors) +
  my_theme +
  coord_equal()

ggsave(TSNE_clusters + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))),
       filename = here(plot_dir, "clusters_HC-29.png"), width = 10)


#### compare w/ mb k-means ####
load(here("processed-data", "03_build_sce","km_res.Rdata"), verbose = TRUE)
k_list <- seq(5, 50) ## keep index
sce$kmeans<- as.factor(paste0("mbk", numform::f_pad_zero(km_res[[which(k_list==29)]]$Clusters)))

cluster_compare <- table(sce$collapsedCluster, sce$kmeans)

library("pheatmap")

cluster_compare2 <- cluster_compare
cluster_compare2[cluster_compare > 5000] <- 5000

png(here(plot_dir, "cluster_compare_heatmap.png"),height = 800, width = 800)
pheatmap(cluster_compare2)
dev.off()


#### Marker Genes ####
# Just for logcounts
sce <- batchelor::multiBatchNorm(sce, batch=sce$Sample)

# load Mathy's markers
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/markers.rda", verbose = TRUE)
# markers.mathys.custom
all(unlist(markers.mathys.custom) %in% rowData(sce)$gene_name)

rownames(sce) <- rowData(sce)$gene_name

source("my_plotExpression.R")

## Mathys
my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.custom,
               cat = "collapsedCluster",
               fill_colors = cluster_colors,
               pdf_fn = here(plot_dir, "HC_mathys_markers.pdf"))

## Tran 
dlpfc_markers <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv")
dlpfc_markers_list <- as.list(dlpfc_markers[1:4,grepl("_1vAll", colnames(dlpfc_markers))])

my_plotMarkers(sce = sce, 
               marker_list = dlpfc_markers_list,
               cat = "collapsedCluster",
               fill_colors = cluster_colors,
               pdf_fn = here(plot_dir, "HC_Tran_markers.pdf"))

## Mean Ratio Top Markers
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

mr_list <- marker_stats %>% 
  ungroup() %>%
  arrange(cellType.target) %>%
  filter(rank_ratio < 5) %>%
  select(rank_ratio, cellType.target, Symbol) %>%
  tidyr::pivot_wider(names_from = 'cellType.target', values_from = 'Symbol') %>%
  select(-rank_ratio) %>%
  as.list()

my_plotMarkers(sce = sce, 
               marker_list = mr_list,
               cat = "collapsedCluster",
               fill_colors = cluster_colors,
               pdf_fn = here(plot_dir, "HC_MeanRatio_markers.pdf"))

my_plotClusterMarkers(sce = sce,
                      marker_list = markers.mathys.custom,
                      cat = "collapsedCluster",
                      pdf_fn = here(plot_dir, "pcm_test2.pdf"))

my_plotClusterMarkers(sce = sce,
                      marker_list = markers.mathys.custom,
                      cat = "kmeans",
                      pdf_fn = here(plot_dir, "pcm_test2.pdf"))

my_plotClusterMarkers(sce = sce,
                      marker_list = mr_list,
                      cat = "collapsedCluster",
                      pdf_fn = here(plot_dir, "pcm_MR.pdf"))

pdf(here(plot_dir, "pcm_test.pdf"), height=6, width=8)
my_plotExpression_flip(sce[,sce$collapsedCluster == "HC01"],
               marker_list = markers.mathys.custom, 
               assay = "logcounts")
dev.off()

## check out annotation 
anno_k <- read.csv(here("processed-data", "03_build_sce", "DLPFC_k29_anno.csv"))

anno_k2 <- anno_k %>% 
  arrange(broad) %>%
  group_by(broad) %>% 
  mutate(ci = as.character(row_number()),
         ct = paste0(broad,'_',ci),
         cellType = ifelse(n()>1,"multi",broad))

anno_k2 %>% select(broad,ci,ct,cellType)

anno_k %>% group_by(broad) %>% mutate(cellType = paste0(broad,"_",rownumber()))
table(anno_k$broad)

anno_hc <- read.csv(here("processed-data", "03_build_sce", "DLPFC_HC_anno.csv"))
table(anno_hc$broad)

sce$cellType_k <- factor(anno_k$broad[match(sce$kmeans, anno_k$cluster)])
table(sce$cellType_k)

sce$cellType_hc <- factor(anno_hc$broad[match(sce$collapsedCluster, anno_hc$cluster)])
table(sce$cellType_hc)
table(sce$collapsedCluster)

ct_tab <- table(sce$cellType_hc,sce$cellType_k)

png(here(plot_dir, "cellType_compare_heatmap.png"),height = 800, width = 800)
pheatmap(ct_tab)
dev.off()

