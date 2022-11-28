library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("dendextend")
library("dynamicTreeCut")
library("here")
library("sessioninfo")
library("dplyr")
library("numform")
library("pheatmap")

# Plotting set up
my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "03_build_sce", "08_cluster_explore")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data", "03_build_sce", "cell_type_colors.Rdata"), verbose = TRUE)

## Load sce object
load(here("processed-data", "03_build_sce", "sce_harmony_Sample.Rdata"), verbose = TRUE)
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
round_prelimClusters <- table(sce$prelimCluster, sce$round)
round_prelimClusters[which(rowSums(round_prelimClusters != 0) == 1), ]
## 270 clusters from round1
length(round_prelimClusters[which(rowSums(round_prelimClusters != 0) == 1), ])

sample_prelimClusters <- table(sce$prelimCluster, sce$Sample)
sample_prelimClusters[which(rowSums(sample_prelimClusters != 0) == 1), ]
# Br6432_ant 2 exclusive, Br6471_ant 7 exclusive

# # rbind the ref.sampleInfo[.rev]
# ref.sampleInfo <- rbind(ref.sampleInfo, ref.sampleInfo.rev)

#### check doublet score for each prelim clust ####
clusIndexes <- splitit(sce$prelimCluster)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii) {
    median(sce$doubletScore[ii])
})

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

## Prelim cluster 156 Over 5 but does not drive clustering
(doublet_clust <- prelimCluster.medianDoublet[prelimCluster.medianDoublet > 5])
# 156
# 8.395549
sum(sce$prelimCluster == 156)
# [1] 68
summary(sce$doubletScore[sce$prelimCluster == 156])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3844  1.6301  8.3955  7.9715 13.7577 18.0338

table(sce$prelimCluster)[names(doublet_clust)]
# 156
# 68

## Cluster in 07_hierachiacal_cluster.R - need to organize
load(here("processed-data", "03_build_sce", "HC_dend.Rdata"), verbose = TRUE)
# dend
# tree.clusCollapsed
# dist.clusCollapsed

## Try different cut heights
h_list <- seq(100, 2000, 25)

h_n_clust <- sapply(h_list, function(h) {
    clust.treeCut <- cutreeDynamic(tree.clusCollapsed,
        distM = as.matrix(dist.clusCollapsed),
        minClusterSize = 2, deepSplit = 1, cutHeight = h
    )
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

clust.treeCut <- cutreeDynamic(tree.clusCollapsed,
    distM = as.matrix(dist.clusCollapsed),
    minClusterSize = 2, deepSplit = 1, cutHeight = chosen_cut_height
)


table(clust.treeCut)
unname(clust.treeCut[order.dendrogram(dend)])
unname(clust.treeCut[order.dendrogram(dend)]) == 0


# Add new labels to those prelimClusters cut off
## just define as a max cluster for now
if (any(clust.treeCut[order.dendrogram(dend)] == 0)) {
    max_clust <- max(clust.treeCut)
    clust.treeCut[order.dendrogram(dend)][which(clust.treeCut[order.dendrogram(dend)] == 0)] <- max_clust + 1

    # 'Re-write', since there are missing numbers
    clust.treeCut[order.dendrogram(dend)] <- as.numeric(as.factor(clust.treeCut[order.dendrogram(dend)]))
    cluster_colors[as.character(max_clust + 1)] <- "black"
}

## Define HC color pallet
names(clust.treeCut) <- paste0("HC", numform::f_pad_zero(names(clust.treeCut)))
cluster_colors <- DeconvoBuddies::create_cell_colors(cell_types = sort(unique(names(clust.treeCut))), pallet = "gg")

labels_colors(dend) <- cluster_colors[clust.treeCut[order.dendrogram(dend)]]

# Print for future reference
pdf(here(plot_dir, "dend_cut.pdf"), height = 11)
par(cex = 0.3, font = 1)
plot(dend, main = "DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = chosen_cut_height, lty = 2)
dev.off()

# png(here(plot_dir, "dend_cut.png"),width = 1000, height = 2000)
png(here(plot_dir, "dend_cut.png"), height = 10, width = 6, unit = "in", res = 400)
par(cex = 0.3, font = 1)
plot(dend, main = "DLPFC prelim-kNN-cluster relationships with collapsed assignments", horiz = TRUE)
abline(v = chosen_cut_height, lty = 2)
dev.off()

# Make reference for new cluster assignment
clusterRefTab.dlpfc <- data.frame(
    origClust = order.dendrogram(dend),
    merged = names(clust.treeCut)[order.dendrogram(dend)]
)
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
tail(table(sce$prelimCluster, sce$collapsedCluster), 20)


#### Plot clusters in TSNE ####
TSNE_clusters <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = collapsedCluster)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cluster_colors) +
    my_theme +
    coord_equal()

ggsave(
    TSNE_clusters +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_HC-29_clusters.png"), width = 10
)


#### compare w/ mb k-means ####
load(here("processed-data", "03_build_sce", "km_res.Rdata"), verbose = TRUE)
k_list <- seq(5, 50) ## keep index
sce$kmeans <- as.factor(paste0("mbk", numform::f_pad_zero(km_res[[which(k_list == 29)]]$Clusters)))
table(sce$kmeans)


#### QC prelim clusters with poor signal downstream ####
# clusters from Oligo_01 32, 36, 180, 229, and 247 - poor signal

problem_clusters <- c(32, 36, 180, 229, 247)
table(sce$prelimCluster %in% problem_clusters)
# FALSE  TRUE
# 56447 21157

## Add QC_remove as a collapse level - signaling it will be dropped
levels(sce$collapsedCluster) <- c(levels(sce$collapsedCluster), "QC_remove")

## label probelm clusters "QC_remove"
sce$collapsedCluster[sce$prelimCluster %in% problem_clusters] <- "QC_remove"
table(sce$collapsedCluster)

## Replot TSNE/UMAP, need to add color for QC_remove
cluster_colors <- c(cluster_colors, c(`QC_remove` = "grey"))

TSNE_clusters <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = collapsedCluster)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cluster_colors) +
    my_theme +
    coord_equal()

ggsave(
    TSNE_clusters +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_HC-29qc_clusters.png"), width = 10
)

####  get logcounts ####
sce <- logNormCounts(sce)
rownames(sce) <- rowData(sce)$gene_name


#### Load data ####
# load(here("processed-data", "sce","sce_DLPFC.Rdata"), verbose = TRUE)


#### Marker Genes ####
# load Mathy's markers
load(here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"), verbose = TRUE)
# markers.mathys.tran
all(unlist(markers.mathys.tran) %in% rowData(sce)$gene_name)
# [1] TRUE

source("my_plotExpression.R")

## Mathys
my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "collapsedCluster",
    fill_colors = cluster_colors,
    pdf_fn = here(plot_dir, "markers_mathys_hc_29.pdf")
)

## iWantHue 29 colors
load(here("processed-data", "03_build_sce", "color_palletes.Rdata"), verbose = TRUE)

my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "kmeans",
    fill_colors = iWantHue_k29,
    pdf_fn = here(plot_dir, "markers_mathys_mb_kmeans_29.pdf")
)

#### Other markers to check ####
# dlpfc_markers <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv")
# dlpfc_markers_list <- as.list(dlpfc_markers[1:4, grepl("_1vAll", colnames(dlpfc_markers))])
#
# my_plotMarkers(
#     sce = sce,
#     marker_list = dlpfc_markers_list,
#     cat = "collapsedCluster",
#     fill_colors = cluster_colors,
#     pdf_fn = here(plot_dir, "HC_Tran_markers.pdf")
# )
#
# ## Mean Ratio Top Markers
# load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)
#
# mr_list <- marker_stats %>%
#     ungroup() %>%
#     arrange(cellType.target) %>%
#     filter(rank_ratio < 5) %>%
#     select(rank_ratio, cellType.target, Symbol) %>%
#     tidyr::pivot_wider(names_from = "cellType.target", values_from = "Symbol") %>%
#     select(-rank_ratio) %>%
#     as.list()
#
# my_plotMarkers(
#     sce = sce,
#     marker_list = mr_list,
#     cat = "collapsedCluster",
#     fill_colors = cluster_colors,
#     pdf_fn = here(plot_dir, "HC_MeanRatio_markers.pdf")
# )


#### Assign annotation ####
## mbkm annotation
anno_k <- read.csv(here("processed-data", "03_build_sce", "DLPFC_k29_anno.csv"))
table(anno_k$broad)
# Astro       drop  EndoMural      Excit      Inhib MicroOligo      Oligo        OPC
# 1          4          1         13          6          1          2          1

## Add index to multiple cell types
anno_k2 <- anno_k %>%
    group_by(broad) %>%
    mutate(
        ct = paste0(broad, "_", f_pad_zero(row_number(), width = 2)),
        n = n()
    ) %>%
    ungroup() %>%
    mutate(cellType = as.factor(ifelse(n > 1, ct, broad)))

anno_k2 %>%
    select(broad, ct, cellType) %>%
    arrange(cellType)

anno_k2 %>% filter(n == 1)

## HC annotation
anno_hc <- read.csv(here("processed-data", "03_build_sce", "DLPFC_HC_anno.csv"))
table(anno_hc$broad)
# Astro      drop EndoMural     Excit     Inhib     Micro     Oligo       OPC
# 1         1         2        15         6         1         3         1


anno_hc2 <- anno_hc %>%
    group_by(broad) %>%
    mutate(
        ct = paste0(broad, "_", f_pad_zero(row_number(), width = 2)),
        n = n()
    ) %>%
    ungroup() %>%
    mutate(cellType = as.factor(ifelse(n > 1, ct, broad)))

order_cell_types <- function(cell_types) {
    neun_mask <- grepl("Excit|Inhib", cell_types)
    neun_ct <- cell_types[neun_mask]
    glia_ct <- cell_types[!neun_mask]
    ordered_cell_types <- c(sort(glia_ct), sort(neun_ct))
    return(ordered_cell_types)
}


## Assign to sce
sce$cellType_broad_k <- factor(anno_k$broad[match(sce$kmeans, anno_k$cluster)])
sce$cellType_k <- factor(anno_k2$cellType[match(sce$kmeans, anno_k2$cluster)])

table(sce$cellType_broad_k)
# Astro       drop  EndoMural      Excit      Inhib MicroOligo      Oligo        OPC
# 3557         23       1330      21233      10413       5541      33716       1791

table(sce$cellType_k)
# Astro    drop_01    drop_02    drop_03    drop_04  EndoMural   Excit_01   Excit_02   Excit_03   Excit_04
# 3557         11          3          2          7       1330       5879       1058        548       1188
# Excit_05   Excit_06   Excit_07   Excit_08   Excit_09   Excit_10   Excit_11   Excit_12   Excit_13   Inhib_01
# 82        633       3198       2851        896       1581       1732       1292        295        132
# Inhib_02   Inhib_03   Inhib_04   Inhib_05   Inhib_06 MicroOligo   Oligo_01   Oligo_02        OPC
# 461       2242       6018       1311        249       5541      28085       5631       1791


hc_levels <- order_cell_types(levels(anno_hc2$cellType))
hc_broad_levels <- order_cell_types(unique(anno_hc2$broad))

sce$cellType_broad_hc <- factor(anno_hc$broad[match(sce$collapsedCluster, anno_hc$cluster)], levels = hc_broad_levels)
sce$cellType_hc <- factor(anno_hc2$cellType[match(sce$collapsedCluster, anno_hc2$cluster)], levels = hc_levels)
table(sce$cellType_broad_hc)
# Astro      drop EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979     21157      2157      1601     10894      1940     24809     11067

table(sce$cellType_hc)
# Astro         drop EndoMural_01 EndoMural_02        Micro     Oligo_01     Oligo_02     Oligo_03          OPC
# 3979        21157          446         1711         1601         1868         4732         4294         1940
# Excit_01     Excit_02     Excit_03     Excit_04     Excit_05     Excit_06     Excit_07     Excit_08     Excit_09
# 7927         2487         1309         2171         2532          329          334         1463         2561
# Excit_10     Excit_11     Excit_12     Excit_13     Excit_14     Excit_15     Inhib_01     Inhib_02     Inhib_03
# 1079          482          420         1567           82           66         5366         1267         1310
# Inhib_04     Inhib_05     Inhib_06
# 565         1192         1367

## Check out prop
(prop_k <- 100 * round(table(sce$cellType_broad_k) / ncol(sce), 3))
# Astro        drop  EndoMural       Excit       Inhib MicroOligo       Oligo         OPC
# 4.6         0.0         1.7        27.4        13.4         7.1        43.4         2.3

(prop_hc <- 100 * round(table(sce$cellType_broad_hc) / ncol(sce), 3))
# Astro      drop EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 5.1      27.3       2.8       2.1      14.0       2.5      32.0      14.3

as.data.frame(prop_k) %>% full_join(as.data.frame(prop_hc) %>% rename(Freq_HC = Freq))
#         Var1 Freq Freq_HC
# 1      Astro  4.6     5.1
# 2       drop  0.0    27.3
# 3  EndoMural  1.7     2.8
# 4      Excit 27.4    32.0
# 5      Inhib 13.4    14.3
# 6 MicroOligo  7.1      NA
# 7      Oligo 43.4    14.0 * QC looses lots of Oligo nuc
# 8        OPC  2.3     2.5
# 9      Micro   NA     2.1

## Tran, Maynard DLPFC dataset for refrence
# Astro      Excit      Inhib Macrophage      Micro      Mural      Oligo        OPC      Tcell
# 7.0       21.3       14.1        0.1        3.5        0.2       48.7        5.1        0.1

## Export all cell types for ref
(ct <- sort(unique(as.character(c(sce$cellType_k, sce$cellType_hc)))))
cat(ct, file = here("processed-data", "03_build_sce", "cell_types.txt"), sep = "\n")

#### Match colData to other DLPFC datasets ####

table(sce$Sample) # good!
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant
# 3101        5911        2861        2723        5205        4282        3898        4067        3059        3212
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid
# 4724        4004        4275        4707        4020        4997        2661        5774        4123

table(sce$region) # need to swap to sentence case and change col name to "Position"
# anterior    middle posterior
# 28716     31974     16914

sce$region <- stringr::str_to_title(sce$region)


table(sce$region_short) # need to change colname to "pos"
# ant   mid  post
# 28716 31974 16914
table(sce$file_id) # need to change colname to "SAMPLE_ID"
# 10c_k  11c_k  12c_k  13c_k  14c_k  15c_k  16c_k  17c_k  18c_k   1c-k   2c-k   3c-k   4c_k   5c_k   6c_k   7c_k   8c_k
# 4997   5205   4067   4123   3898   4282   4707   4020   5774   3101   3059   3212   4004   4275   4724   2661   2861
# 9c_k round0
# 5911   2723

## Fix colnames
cn <- names(colData(sce))
cn[match(c("subject", "region", "region_short", "file_id"), cn)] <- c("BrNum", "Position", "pos", "SAMPLE_ID")
colnames(colData(sce)) <- cn

sce$Position <- tools::toTitleCase(sce$Position)

table(sce$BrNum, sce$Position)
# Anterior Middle Posterior
# Br2720        0   3101      5911
# Br2743     2861   2723         0
# Br3942     5205   4282         0
# Br6423     3898      0      4067
# Br6432     3059      0         0
# Br6471     3212   4724         0
# Br6522        0   4004      4275
# Br8325     4707   4020         0
# Br8492        0   4997      2661
# Br8667     5774   4123         0


#### DROP QC'ed HC cells ####

# sce <- sce[, sce$cellType_hc != "drop"]
#
# sce$cellType_hc <- droplevels(sce$cellType_hc)
# levels(sce$cellType_hc)

#### SAVE ANNOTATED SCE ####
## add cell type colors as metadata
metadata(sce)$cell_type_colors <- cell_type_colors
metadata(sce)$cell_type_colors_broad <- cell_type_colors_broad

save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('08_cluster_annotation', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript 08_cluster_annotation.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
