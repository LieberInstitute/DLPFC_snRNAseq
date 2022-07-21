library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")
library("dplyr")
library("numform")

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

## Cluster in 07_hierachiacal_cluster.R - need to organize 
load(here("processed-data", "03_build_sce", "HC_dend.Rdata"), verbose = TRUE)

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
pdf(here(plot_dir, "dend_cut.pdf"), height = 11)
par(cex=0.3, font=1)
plot(dend, main="DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = chosen_cut_height, lty = 2)
dev.off()

# png(here(plot_dir, "dend_cut.png"),width = 1000, height = 2000)
png(here(plot_dir, "dend_cut.png"), height = 10, width = 6, unit = "in", res = 400)
par(cex=0.3, font=1)
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
# HC01  HC02  HC03  HC04  HC05  HC06  HC07  HC08  HC09  HC10  HC11  HC12  HC13  HC14  HC15  HC16  HC17  HC18  HC19  HC20 
# 23025  7927  2487  5366  1309  2171  4732  1267  2532   329  4294  1940   334  1463  1310   565  2561  3979  1192  1367 
# HC21  HC22  HC23  HC24  HC25  HC26  HC27  HC28  HC29 
# 1079   482   420  1601  1567   446  1711    82    66 

n_clusters <- length(levels(sce$collapsedCluster))
# Print some visualizations:
tail(table(sce$prelimCluster, sce$collapsedCluster),20)


#### Plot clusters in TSNE ####
TSNE_clusters <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=collapsedCluster)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cluster_colors) +
  my_theme +
  coord_equal()

ggsave(TSNE_clusters + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))),
       filename = here(plot_dir, "TSNE_HC-29_clusters.png"), width = 10)


# get logcounts
sce <- logNormCounts(sce)
rownames(sce) <- rowData(sce)$gene_name

#### compare w/ mb k-means ####
load(here("processed-data", "03_build_sce","km_res.Rdata"), verbose = TRUE)
k_list <- seq(5, 50) ## keep index
sce$kmeans<- as.factor(paste0("mbk", numform::f_pad_zero(km_res[[which(k_list==29)]]$Clusters)))
table(sce$kmeans)


#### Load data ####
load(here("processed-data", "03_build_sce","cell_type_colors.Rdata"), verbose = TRUE)
# cell_type_colors_broad
# cell_type_colors
# load(here("processed-data", "sce","sce_DLPFC.Rdata"), verbose = TRUE)


#### Marker Genes ####
# load Mathy's markers
load(here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"), verbose = TRUE)
# markers.mathys.tran
all(unlist(markers.mathys.tran) %in% rowData(sce)$gene_name)
# [1] TRUE

source("my_plotExpression.R")

## Mathys
my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "collapsedCluster",
               fill_colors = cluster_colors,
               pdf_fn = here(plot_dir, "HC_mathys_markers.pdf"))

my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "kmeans",
               fill_colors = iWantHue_k29,
               pdf_fn = here(plot_dir, "mb_kmeans_29_Tran_markers.pdf"))

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
                      marker_list = mr_list,
                      cat = "collapsedCluster",
                      pdf_fn = here(plot_dir, "pcm_MR.pdf"))

pdf(here(plot_dir, "pcm_test.pdf"), height=6, width=8)
my_plotExpression_flip(sce[,sce$kmeans == "mbk05"],
               marker_list = markers.mathys.tran, 
               assay = "logcounts")
dev.off()

#### Assign annotation ####
## mbkm annotation
anno_k <- read.csv(here("processed-data", "03_build_sce", "DLPFC_k29_anno.csv"))
table(anno_k$broad)
# Astro        drop  Endo.Mural       Excit       Inhib Micro.Oligo       Oligo         OPC 
# 1           4           1          13           6           1           2           1 

## Add index to multiple cell types
anno_k2 <- anno_k %>% 
  group_by(broad) %>% 
  mutate(ct = paste0(broad,'_',f_pad_zero(row_number(), width = 2)),
         n = n()) %>%
  ungroup() %>%
  mutate(cellType = as.factor(ifelse(n>1,ct,broad)))

anno_k2 %>% select(broad,ct,cellType) %>% arrange(cellType)
anno_k2 %>% filter(n == 1)

## HC annotation 
anno_hc <- read.csv(here("processed-data", "03_build_sce", "DLPFC_HC_anno.csv"))
table(anno_hc$broad)
# Astro Endo.Mural      Excit      Inhib      Micro      Oligo        OPC 
# 1          2         15          6          1          3          1 


anno_hc2 <- anno_hc %>% 
  group_by(broad) %>% 
  mutate(ct = paste0(broad,'_',f_pad_zero(row_number(), width = 2)),
         n = n()) %>%
  ungroup() %>%
  mutate(cellType = as.factor(ifelse(n>1,ct,broad)))


## Assign to sce
sce$cellType_broad_k <- factor(anno_k$broad[match(sce$kmeans, anno_k$cluster)])
sce$cellType_k <- factor(anno_k2$cellType[match(sce$kmeans, anno_k2$cluster)])

table(sce$cellType_broad_k)
# Astro        drop  Endo.Mural       Excit       Inhib Micro.Oligo       Oligo         OPC 
# 3557          23        1330       21233       10413        5541       33716        1791

table(sce$cellType_k)
# Astro     drop_01     drop_02     drop_03     drop_04  Endo.Mural    Excit_01    Excit_02    Excit_03    Excit_04 
# 3557          11           3           2           7        1330        5879        1058         548        1188 
# Excit_05    Excit_06    Excit_07    Excit_08    Excit_09    Excit_10    Excit_11    Excit_12    Excit_13    Inhib_01 
# 82         633        3198        2851         896        1581        1732        1292         295         132 
# Inhib_02    Inhib_03    Inhib_04    Inhib_05    Inhib_06 Micro.Oligo    Oligo_01    Oligo_02         OPC 
# 461        2242        6018        1311         249        5541       28085        5631        1791 


sce$cellType_broad_hc <- factor(anno_hc$broad[match(sce$collapsedCluster, anno_hc$cluster)])
sce$cellType_hc <- factor(anno_hc2$cellType[match(sce$collapsedCluster, anno_hc2$cluster)])
table(sce$cellType_broad_hc)
# Astro Endo.Mural      Excit      Inhib      Micro      Oligo        OPC 
# 3979       2157      24809      11067       1601      32051       1940 

table(sce$cellType_hc)
# Endo.Mural_01 Endo.Mural_02      Excit_01      Excit_02      Excit_03      Excit_04      Excit_05      Excit_06 
#           446          1711          7927          2487          1309          2171          2532           329 
# Excit_07      Excit_08      Excit_09      Excit_10      Excit_11      Excit_12      Excit_13      Excit_14 
#      334          1463          2561          3979          1079           482           420          1567 
# Excit_15      Excit_16      Inhib_01      Inhib_02      Inhib_03      Inhib_04      Inhib_05      Inhib_06 
#       82            66          5366          1267          1310           565          1192          1367 
# Micro      Oligo_01      Oligo_02      Oligo_03           OPC 
#  1601         23025          4732          4294          1940 

## Check out prop
(prop_k <- 100*round(table(sce$cellType_broad_k)/ncol(sce),3))
# Astro        drop  Endo.Mural       Excit       Inhib Micro.Oligo       Oligo         OPC 
# 4.6         0.0         1.7        27.4        13.4         7.1        43.4         2.3 

(prop_hc <- 100*round(table(sce$cellType_broad_hc)/ncol(sce),3))
# Astro Endo.Mural      Excit      Inhib      Micro      Oligo        OPC 
# 5.1        2.8       32.0       14.3        2.1       41.3        2.5 

as.data.frame(prop_k) %>% full_join(as.data.frame(prop_hc) %>% rename(Freq_HC = Freq))
#          Var1 Freq Freq_HC
# 1       Astro  4.6     5.1
# 2        drop  0.0      NA
# 3  Endo.Mural  1.7     2.8
# 4       Excit 27.4    32.0
# 5       Inhib 13.4    14.3
# 6 Micro.Oligo  7.1      NA
# 7       Oligo 43.4    41.3
# 8         OPC  2.3     2.5
# 9       Micro   NA     2.1

## Tran, Maynard DLPFC dataset for refrence
# Astro      Excit      Inhib Macrophage      Micro      Mural      Oligo        OPC      Tcell 
# 7.0       21.3       14.1        0.1        3.5        0.2       48.7        5.1        0.1 

save(sce, file = here("processed-data", "sce","sce_DLPFC.Rdata"))
# load(here("processed-data", "sce","sce_DLPFC.Rdata"), verbose = TRUE)

#### Replot with Annotations ####
(ct <- sort(unique(as.character(c(sce$cellType_k,sce$cellType_hc)))))
cat(ct, file = here("processed-data","03_build_sce","cell_types.txt"), sep = "\n")
all(ct %in% names(cell_type_colors)) ## if not revisit cell_colors

## Plot clusters in TSNE
TSNE_HC_cellTypes <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=cellType_hc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
  my_theme +
  coord_equal()

ggsave(TSNE_HC_cellTypes + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))),
       filename = here(plot_dir, "TSNE_HC-29_cellType.png"), width = 10)


TSNE_HC_cellTypes <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2)) +
  geom_point(aes(colour=cellType_broad_hc), size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors[levels(sce$cellType_broad_hc)], drop = TRUE) +
  my_theme +
  coord_equal()

ggsave(TSNE_HC_cellTypes + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))),
       filename = here(plot_dir, "TSNE_HC_broad_cellType.png"), width = 10)



TSNE_km_cellTypes <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=cellType_k)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
  my_theme +
  coord_equal()

ggsave(TSNE_km_cellTypes + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))),
       filename = here(plot_dir, "TSNE_mbkm-29_cellType.png"), width = 10)


## markers
my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "HC_mathys_markers_ct.pdf"))

my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "cellType_k",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "mb_kmeans_29_mathys_markers_ct.pdf"))


#### Compare annotaitons ####
library("pheatmap")

table(sce$kmeans)
table(sce$collapsedCluster)

cluster_compare <- table(sce$kmeans, sce$collapsedCluster)
cluster_compare_prop <- sweep(cluster_compare,2,colSums(cluster_compare),`/`)

cluster_compare > 0

hc_counts <- table(sce$kmeans)

hc_anno <- anno_hc %>% select(cluster, broad) %>% tibble::column_to_rownames("cluster") %>% mutate(n = table(sce$collapsedCluster))
km_anno <- anno_k %>% select(cluster, broad) %>% tibble::column_to_rownames("cluster") %>% mutate(n = table(sce$kmeans))

broad_pallet <- cell_type_colors_broad[unique(c(hc_anno$broad, km_anno$broad))] 

png(here(plot_dir, "cluster_compare_heatmap.png"),height = 800, width = 800)
pheatmap(cluster_compare_prop,
         annotation_col= hc_anno,
         annotation_row = km_anno,
         annotation_colors = list(broad = temp_pallet)
         )
dev.off()

## Jaccard indicies
library(bluster)
jacc.mat <- linkClustersMatrix(sce$kmeans, sce$collapsedCluster)

png(here(plot_dir, "cluster_compare_heatmap_jacc.png"),height = 800, width = 800)
pheatmap(jacc.mat,
         annotation_col= hc_anno,
         annotation_row = km_anno,
         annotation_colors = list(broad = temp_pallet)
)
dev.off()

## best cluster pairing
best <- max.col(jacc.mat, ties.method="first")
DataFrame(
  Cluster=rownames(jacc.mat), 
  Corresponding=colnames(jacc.mat)[best], 
  Index=jacc.mat[cbind(seq_len(nrow(jacc.mat)), best)]
)

## compare cell types
ct_tab <- table(sce$cellType_k, sce$cellType_hc)
ct_tab_prop <- sweep(ct_tab,2,colSums(ct_tab),`/`)

hc_ct_anno <- table(sce$cellType_hc) %>% as.data.frame() %>% rename(ct = Var1, n = Freq)
rownames(hc_ct_anno) <- hc_ct_anno$ct

km_ct_anno <- table(sce$cellType_k) %>% as.data.frame() %>% rename(ct = Var1, n = Freq)
rownames(km_ct_anno) <- km_ct_anno$ct

png(here(plot_dir, "cellType_compare_heatmap.png"),height = 800, width = 800)
pheatmap(ct_tab_prop,
         annotation_col = hc_ct_anno,
         annotation_row = km_ct_anno,
         annotation_colors = list(ct = cell_type_colors)
         )
dev.off()

jacc.mat.ct <- linkClustersMatrix(sce$cellType_k, sce$cellType_hc)

png(here(plot_dir, "cellType_compare_heatmap_jacc.png"),height = 800, width = 800)
pheatmap(jacc.mat.ct,
         annotation_col = hc_ct_anno,
         annotation_row = km_ct_anno,
         annotation_colors = list(ct = cell_type_colors)
)
dev.off()

## compare broad cell types

(ctb_tab <- table(sce$cellType_broad_k, sce$cellType_broad_hc))
ctb_tab_prop <- sweep(ctb_tab,2,colSums(ctb_tab),`/`)

hc_broad_anno <- table(sce$cellType_broad_hc) %>% as.data.frame() %>% rename(broad = Var1, n = Freq)
rownames(hc_broad_anno) <- hc_broad_anno$broad

km_broad_anno <- table(sce$cellType_broad_k) %>% as.data.frame() %>% rename(broad = Var1, n = Freq)
rownames(km_broad_anno) <- km_broad_anno $broad
  
png(here(plot_dir, "cellType_broad_compare_heatmap.png"),height = 800, width = 800)
pheatmap(ctb_tab_prop,
         annotation_col = hc_broad_anno,
         annotation_row = km_broad_anno,
         annotation_colors = list(broad = broad_pallet))
dev.off()

jacc.mat.broad <- linkClustersMatrix(sce$cellType_broad_k, sce$cellType_broad_hc)

png(here(plot_dir, "cellType_broad_compare_heatmap_jacc.png"),height = 800, width = 800)
pheatmap(jacc.mat.broad ,
         annotation_col = hc_broad_anno,
         annotation_row = km_broad_anno,
         annotation_colors = list(broad = broad_pallet))
dev.off()

## Adjusted Rand Index
# 0.5 corresponds to “good” similarity
pairwiseRand(sce$kmeans, sce$collapsedCluster, mode="index")
# [1] 0.5875801
pairwiseRand(sce$cellType_broad_k, sce$cellType_broad_hc, mode="index")
# [1] 0.5791338

#### Explore UMI ####
UMI_ct_k <- ggcells(sce, mapping=aes(x=cellType_k, y=sum, fill=cellType_k)) +
  geom_boxplot()+
  scale_fill_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
  my_theme +
  theme(legend.position = "None",axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1))

ggsave(UMI_ct_k,filename = here(plot_dir, "UMI_mbkm-29_cellType.png"), width = 10)


UMI_ct_hc <- ggcells(sce, mapping=aes(x=cellType_hc, y=sum, fill=cellType_hc)) +
  geom_boxplot()+
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(legend.position = "None",axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1))

ggsave(UMI_ct_hc,filename = here(plot_dir, "UMI_HC_cellType.png"), width = 10)


#### Explore doublet scores
dbs_ct_k <- ggcells(sce, mapping=aes(x=cellType_k, y=doubletScore, fill=cellType_k)) +
  geom_boxplot()+
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(legend.position = "None",axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1))

ggsave(dbs_ct_k,filename = here(plot_dir, "doubletScore_mbkm-29_cellType.png"), width = 10)


dbs_ct_hc <- ggcells(sce, mapping=aes(x=cellType_hc, y=doubletScore, fill=cellType_hc)) +
  geom_boxplot()+
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(legend.position = "None",axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1))

ggsave(dbs_ct_hc,filename = here(plot_dir, "doubletScore_HC_cellType.png"), width = 10)


cellType.idx <- splitit(sce$cellType_hc)
#sapply(c("Excit", "Inhib", "MSN"), function(x){grep(x, names(cellType.idx))})

sapply(cellType.idx, function(x){quantile(sce$doubletScore[x])})[ ,order(sapply(cellType.idx, function(x){quantile(sce$doubletScore[x])["50%"]}))]
#       Excit_15    Micro  Oligo_03 Excit_09     Astro       OPC Excit_10   Oligo_02 Excit_13  Excit_11 Endo.Mural_02
# 0%   0.0086180 0.000000  0.000000  0.00000  0.000000  0.000000 0.000000  0.0000000 0.000000  0.050708      0.000000
# 25%  0.0348360 0.032112  0.042816  0.04134  0.057876  0.068640 0.077310  0.0811380 0.103416  0.144880      0.085632
# 50%  0.0751140 0.099216  0.137040  0.13744  0.143964  0.166296 0.180978  0.1896200 0.197182  0.282516      0.300958
# 75%  0.1874415 0.388596  0.568988  0.45220  0.345052  0.639158 0.426949  0.4642535 0.426949  0.476293      0.765917
# 100% 4.4314060 6.411696 15.852624  6.83060 12.939216 10.329944 7.175700 16.9815360 8.621844 24.470232      8.387922
#       Inhib_03  Excit_12 Inhib_02 Endo.Mural_01 Excit_03 Inhib_06 Inhib_05 Excit_05 Excit_04  Inhib_01 Oligo_01
# 0%    0.006864  0.074860 0.035706      0.011902 0.011902 0.074928 0.047608 0.000000 0.035706  0.000000  0.00000
# 25%   0.239052  0.202122 0.255210      0.150456 0.253524 0.250351 0.279440 0.249900 0.368962  0.354956  0.20670
# 50%   0.383932  0.396758 0.419160      0.422282 0.486492 0.489020 0.512789 0.524496 0.608780  0.688180  0.70438
# 75%   0.738796  0.572276 1.123520      0.744120 1.196036 1.137776 1.246423 1.035300 1.277320  1.187472  1.72550
# 100% 12.329442 17.515992 7.654660      4.900264 9.451632 9.283560 9.712032 9.719232 8.319498 31.475180 12.06863
#       Excit_01 Inhib_04 Excit_08 Excit_14 Excit_02 Excit_06 Excit_07
# 0%    0.023800 0.143892 0.019960 0.081284 0.016268 0.199600 0.137636
# 25%   0.476444 0.345158 0.321237 0.477489 0.469060 0.356304 1.230960
# 50%   0.742440 0.759984 0.778440 0.864066 0.902736 1.057420 1.477152
# 75%   1.166396 1.447852 1.533816 1.438500 1.387359 1.712640 1.785096
# 100% 13.872384 9.825768 8.841504 6.457070 8.640852 5.808176 4.430430

sapply(cellType.idx, function(x){quantile(sce$sum[x])})[ ,order(sapply(cellType.idx, function(x){quantile(sce$sum[x])["50%"]}))]
#      Excit_15 Oligo_01 Excit_09 Oligo_03 Micro Endo.Mural_02   Astro Excit_13 Oligo_02 Excit_10 Excit_05      OPC
# 0%        678      220      462      359   350         407.0   276.0    648.0    968.0      533   292.00   397.00
# 25%      1101     1232     1389     1464  1756        2020.5  1474.0   2447.5   3750.0     3298  3077.25  3933.25
# 50%      1218     1907     2121     2211  2633        3140.0  3280.0   3445.0   5113.5     5489  5668.50  5951.00
# 75%      1348     3147     3225     3410  4300        4546.0  7305.5   4680.5   7314.5     8100 10000.00  8510.00
# 100%    15127    27828    50703    36108 15576       32302.0 61842.0  16113.0  49758.0    95913 65772.00 83867.00
#      Excit_14 Endo.Mural_01 Inhib_03 Inhib_06  Inhib_05 Inhib_01 Inhib_04 Excit_11 Excit_12 Excit_03 Inhib_02  Excit_07
# 0%    2853.00       2040.00   641.00   1598.0   1402.00   317.00     1571  4585.00  1215.00     1651   4438.0   2691.00
# 25%   5537.75       6102.75  7002.00   8979.5  10839.25 10929.00    12837 14424.25 13997.75    14598  16177.5  17929.75
# 50%   7893.00       8219.50 12173.00  13173.0  15142.50 16488.00    18805 19989.50 20810.00    20879  23964.0  25357.50
# 75%  17635.75      11667.25 18893.25  20914.5  21695.25 23337.75    26787 28125.75 28306.25    30108  35836.5  39349.00
# 100% 74260.00      44112.00 61063.00  99561.0 135464.00 97656.00    75908 67166.00 59560.00   103322 118342.0 105390.00
#      Excit_04 Excit_01 Excit_08 Excit_06 Excit_02
# 0%     1449.0    327.0    922.0     1437    629.0
# 25%   21769.5  22253.5  19554.5    24773  27708.5
# 50%   31522.0  33633.0  33928.0    37075  42935.0
# 75%   47549.5  51601.5  53335.0    51649  65244.0
# 100% 179194.0 203794.0 188090.0   152333 296099.0

#### Explore layer markers ####
load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_spatial_registration/t_cor_plot_top_genes_k7.rda", verbose = TRUE)
# dat_small
table(dat_small$Layer)
# Layer1 Layer2 Layer3 Layer4 Layer5 Layer6     WM 
# 91     87     65     80     78     83    100 

dat_symb <- rownames(sce)[match(rownames(dat_small),rowData(sce)$gene_id)]
any(is.na(dat_symb))

layer_idx <- splitit(dat_small$Layer)
layer_top4 <- purrr::map(layer_idx, ~dat_symb[.x][1:4])
# $Layer1
# [1] "MT1G"  "VIM"   "FABP7" "MT1F" 
# 
# $Layer2
# [1] "DACT1"   "STXBP6"  "SIPA1L1" "MAN1A1" 
# 
# $Layer3
# [1] "CARTPT"    "BAIAP3"    "ADCYAP1"   "LINC01007"
# 
# $Layer4
# [1] "VAMP1" "NEFH"  "SCN1B" "NGB"  
# 
# $Layer5
# [1] "PCP4"    "CAMK2D"  "SMYD2"   "TRABD2A"
# 
# $Layer6
# [1] "ISLR"  "NR4A2" "DACH1" "NTNG2"
# 
# $WM
# [1] "NDRG1"  "PTP4A2" "AQP1"   "PAQR6"


## markers
my_plotMarkers(sce = sce, 
               marker_list = layer_top4,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "HC_layer_markers_ct.pdf"))

my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "cellType_k",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "mb_kmeans_29_mathys_markers_ct.pdf"))

## Inhib sub-types
markers.Inhib_subtypes = list(
  'inhib_markers pg1' = c("HTR3A", "VIP", "CCK", "NPY", "CRHBP", "CALB2", "PNOC"),
  'inhib_Zhaung2022' = c('SP8', #VIP
                       'KLF5',#SST
                       'LGI2',#PVALB
                       'LAMP5')# LAMP5
)

my_plotMarkers(sce = sce[,grep("Inhib",sce$cellType_hc)], 
               marker_list = markers.Inhib_subtypes,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "HC_inhib_markers_ct.pdf"))

sce_inhib <- sce[,grep("Inhib",sce$cellType_hc)]
sce_inhib$cellType_hc <- droplevels(sce_inhib$cellType_hc)
sce_inhib$collapsedCluster <- droplevels(sce_inhib$collapsedCluster)
table(sce_inhib$cellType_hc, sce_inhib$collapsedCluster)
