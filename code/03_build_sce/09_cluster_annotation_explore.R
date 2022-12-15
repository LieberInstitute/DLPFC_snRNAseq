library("SingleCellExperiment")
library("jaffelab")
library("ggplot2")
library("here")
library("sessioninfo")
library("dplyr")
# library("numform")
library("pheatmap")
library("scater")
library("bluster")
library("patchwork")

# Plotting set up
my_theme <- theme_bw() +
  theme(text = element_text(size = 15))

plot_dir <- here("plots", "03_build_sce", "09_cluster_annotation_explore")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Load Data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

cell_type_colors <- metadata(sce)$cell_type_colors
cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad

#### Replot with Annotations ####
## HC reduce dims will be plotted in 05_explore_sce/01_reduced_dim_plots
# ## Plot clusters in TSNE
# TSNE_HC_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
#   geom_point(size = 0.2, alpha = 0.3) +
#   scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
#   my_theme +
#   coord_equal()
# 
# ggsave(
#   TSNE_HC_cellTypes +
#     guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
#   filename = here(plot_dir, "TSNE_HC-29_cellType.png"), width = 10
# )
# 
# 
# TSNE_HC_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2)) +
#   geom_point(aes(colour = cellType_broad_hc), size = 0.2, alpha = 0.3) +
#   scale_color_manual(values = cell_type_colors_broad[levels(sce$cellType_broad_hc)], drop = TRUE) +
#   my_theme +
#   coord_equal()
# 
# ggsave(
#   TSNE_HC_cellTypes +
#     guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
#   filename = here(plot_dir, "TSNE_HC_broad_cellType.png"), width = 10
# )
# 
# 
# ## UMAP
# UMAP_HC_cellTypes <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType_hc)) +
#   geom_point(size = 0.2, alpha = 0.3) +
#   scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
#   my_theme +
#   coord_equal()
# 
# ggsave(
#   UMAP_HC_cellTypes +
#     guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
#   filename = here(plot_dir, "UMAP_HC_cellType.png"), width = 10
# )
# 
# UMAP_HC_cellTypes <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2)) +
#   geom_point(aes(colour = cellType_broad_hc), size = 0.2, alpha = 0.3) +
#   scale_color_manual(values = cell_type_colors_broad[levels(sce$cellType_broad_hc)], drop = TRUE) +
#   my_theme +
#   coord_equal()
# 
# ggsave(
#   UMAP_HC_cellTypes +
#     guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
#   filename = here(plot_dir, "UMAP_HC_broad_cellType.png"), width = 10
# )

#### Plot KM TSNE + UMAP ####
## jsut plot mbkm here
TSNE_km_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_k)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
  my_theme +
  coord_equal()

ggsave(
  TSNE_km_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
  filename = here(plot_dir, "TSNE_mbkm-29_cellType.png"), width = 10
)


#### Plot markers ####
source("custom_plotExpression.R")
source("my_plotMarkers.R")
load(here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"), verbose = TRUE)

my_plotMarkers(
  sce = sce,
  marker_list = markers.mathys.tran,
  cat = "cellType_hc",
  fill_colors = cell_type_colors,
  pdf_fn = here(plot_dir, "markers_mathys_hc_29_ct.pdf")
)

my_plotMarkers(
  sce = sce,
  marker_list = markers.mathys.tran,
  cat = "cellType_k",
  fill_colors = cell_type_colors,
  pdf_fn = here(plot_dir, "markers_mathys_mb_kmeans_29.pdf")
)


#### Compare annotations ####
table(sce$kmeans)
table(sce$collapsedCluster)

cluster_compare <- table(sce$kmeans, sce$collapsedCluster)
cluster_compare_prop <- sweep(cluster_compare, 2, colSums(cluster_compare), `/`)

hc_anno <- colData(sce) |>
  as.data.frame() |>
  count(collapsedCluster, broad = cellType_broad_hc) |>
  tibble::column_to_rownames("collapsedCluster")

km_anno <- colData(sce) |>
  as.data.frame() |>
  count(kmeans, broad = cellType_broad_k) |>
  tibble::column_to_rownames("kmeans")

broad_pallet <- cell_type_colors_broad[unique(c(hc_anno$broad, km_anno$broad))]

png(here(plot_dir, "heatmap_cluster_compare.png"), height = 800, width = 800)
pheatmap(cluster_compare_prop,
         annotation_col = hc_anno,
         annotation_row = km_anno,
         annotation_colors = list(broad = broad_pallet)
)
dev.off()

## Jaccard indicies
jacc.mat <- linkClustersMatrix(sce$kmeans, sce$collapsedCluster)

png(here(plot_dir, "heatmap_cluster_compare_jacc.png"), height = 800, width = 800)
pheatmap(jacc.mat,
         annotation_col = hc_anno,
         annotation_row = km_anno,
         annotation_colors = list(broad = broad_pallet)
)
dev.off()

## best cluster pairing
best <- max.col(jacc.mat, ties.method = "first")
DataFrame(
  Cluster = rownames(jacc.mat),
  Corresponding = colnames(jacc.mat)[best],
  Index = jacc.mat[cbind(seq_len(nrow(jacc.mat)), best)]
)

## compare cell types
ct_tab <- table(sce$cellType_k, sce$cellType_hc)
ct_tab_prop <- sweep(ct_tab, 2, colSums(ct_tab), `/`)

hc_ct_anno <- table(sce$cellType_hc) |>
  as.data.frame() |>
  rename(ct = Var1, n = Freq)
rownames(hc_ct_anno) <- hc_ct_anno$ct

km_ct_anno <- table(sce$cellType_k) |>
  as.data.frame() |>
  rename(ct = Var1, n = Freq)
rownames(km_ct_anno) <- km_ct_anno$ct

png(here(plot_dir, "heatmap_cellType_compare.png"), height = 800, width = 800)
pheatmap(ct_tab_prop,
         annotation_col = hc_ct_anno,
         annotation_row = km_ct_anno,
         annotation_colors = list(ct = cell_type_colors)
)
dev.off()

jacc.mat.ct <- linkClustersMatrix(sce$cellType_k, sce$cellType_hc)

png(here(plot_dir, "heatmap_cellType_compare_jacc.png"), height = 800, width = 800)
pheatmap(jacc.mat.ct,
         annotation_col = hc_ct_anno,
         annotation_row = km_ct_anno,
         annotation_colors = list(ct = cell_type_colors)
)
dev.off()

## compare broad cell types

(ctb_tab <- table(sce$cellType_broad_k, sce$cellType_broad_hc))
ctb_tab_prop <- sweep(ctb_tab, 2, colSums(ctb_tab), `/`)

hc_broad_anno <- table(sce$cellType_broad_hc) |>
  as.data.frame() |>
  rename(broad = Var1, n = Freq)
rownames(hc_broad_anno) <- hc_broad_anno$broad

km_broad_anno <- table(sce$cellType_broad_k) |>
  as.data.frame() |>
  rename(broad = Var1, n = Freq)
rownames(km_broad_anno) <- km_broad_anno$broad

png(here(plot_dir, "heatmap_cellType_broad_compare.png"), height = 800, width = 800)
pheatmap(ctb_tab_prop,
         annotation_col = hc_broad_anno,
         annotation_row = km_broad_anno,
         annotation_colors = list(broad = broad_pallet)
)
dev.off()

jacc.mat.broad <- linkClustersMatrix(sce$cellType_broad_k, sce$cellType_broad_hc)

png(here(plot_dir, "heatmap_cellType_broad_compare_jacc.png"), height = 800, width = 800)
pheatmap(jacc.mat.broad,
         annotation_col = hc_broad_anno,
         annotation_row = km_broad_anno,
         annotation_colors = list(broad = broad_pallet)
)
dev.off()

## Adjusted Rand Index
# 0.5 corresponds to “good” similarity
pairwiseRand(sce$kmeans, sce$collapsedCluster, mode = "index")
# [1] 0.5875801
pairwiseRand(sce$cellType_broad_k, sce$cellType_broad_hc, mode = "index")
# [1] 0.5791338

#### Explore UMI ####
UMI_ct_k <- ggcells(sce, mapping = aes(x = cellType_k, y = sum, fill = cellType_k)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(UMI_ct_k, filename = here(plot_dir, "qc_UMI_mbkm-29_cellType.png"), width = 10)


UMI_ct_hc <- ggcells(sce, mapping = aes(
  x = reorder(cellType_hc, sum, FUN = median),
  y = sum, fill = cellType_hc
)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Total UMIs")

ggsave(UMI_ct_hc, filename = here(plot_dir, "qc_UMI_HC_cellType.png"), height = 3)
ggsave(UMI_ct_hc, filename = here(plot_dir, "qc_UMI_HC_cellType.pdf"), height = 3)


#### Explore doublet scores
dbs_ct_k <- ggcells(sce, mapping = aes(x = cellType_k, y = doubletScore, fill = cellType_k)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(dbs_ct_k, filename = here(plot_dir, "qc_doubletScore_mbkm-29_cellType.png"), width = 10)


dbs_ct_hc <- ggcells(sce,
                     mapping = aes(
                       x = reorder(cellType_hc, doubletScore, FUN = median),
                       y = doubletScore,
                       fill = cellType_hc
                     )
) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Doublet Score")

ggsave(dbs_ct_hc, filename = here(plot_dir, "qc_doubletScore_HC_cellType.png"), height = 3)
ggsave(dbs_ct_hc, filename = here(plot_dir, "qc_doubletScore_HC_cellType.pdf"), height = 3)

#### Explore Mito rate ####
mito_ct_hc <- ggcells(sce, mapping = aes(
  x = reorder(cellType_hc, subsets_Mito_percent, FUN = median),
  y = subsets_Mito_percent, fill = cellType_hc
)) +
  geom_boxplot() +
  scale_fill_manual(values = cell_type_colors) +
  my_theme +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Percent Mito")

ggsave(mito_ct_hc, filename = here(plot_dir, "qc_mitoRate_HC_cellType.png"), height = 3)
ggsave(mito_ct_hc, filename = here(plot_dir, "qc_mitoRate_HC_cellType.pdf"), height = 3)

## ALL QC plots
ggsave((UMI_ct_hc / dbs_ct_hc / mito_ct_hc), filename = here(plot_dir, "qc_ALL_HC_cellType.png"), height = 7.5)

#### Print out qc metrics ####
cellType.idx <- splitit(sce$cellType_hc)
# sapply(c("Excit", "Inhib", "MSN"), function(x){grep(x, names(cellType.idx))})

sapply(cellType.idx, function(x) {
  quantile(sce$doubletScore[x])
})[, order(sapply(cellType.idx, function(x) {
  quantile(sce$doubletScore[x])["50%"]
}))]
#       Excit_15    Micro  Oligo_03 Excit_09     Astro       OPC Excit_10   Oligo_02 Excit_13  Excit_11 EndoMural_02
# 0%   0.0086180 0.000000  0.000000  0.00000  0.000000  0.000000 0.000000  0.0000000 0.000000  0.050708     0.000000
# 25%  0.0348360 0.032112  0.042816  0.04134  0.057876  0.068640 0.077310  0.0811380 0.103416  0.144880     0.085632
# 50%  0.0751140 0.099216  0.137040  0.13744  0.143964  0.166296 0.180978  0.1896200 0.197182  0.282516     0.300958
# 75%  0.1874415 0.388596  0.568988  0.45220  0.345052  0.639158 0.426949  0.4642535 0.426949  0.476293     0.765917
# 100% 4.4314060 6.411696 15.852624  6.83060 12.939216 10.329944 7.175700 16.9815360 8.621844 24.470232     8.387922
#       Inhib_03  Excit_12 Inhib_02 EndoMural_01 Oligo_01 Excit_03 Inhib_06 Inhib_05 Excit_05 Excit_04  Inhib_01
# 0%    0.006864  0.074860 0.035706     0.011902 0.000000 0.011902 0.074928 0.047608 0.000000 0.035706  0.000000
# 25%   0.239052  0.202122 0.255210     0.150456 0.127732 0.253524 0.250351 0.279440 0.249900 0.368962  0.354956
# 50%   0.383932  0.396758 0.419160     0.422282 0.455001 0.486492 0.489020 0.512789 0.524496 0.608780  0.688180
# 75%   0.738796  0.572276 1.123520     0.744120 1.195782 1.196036 1.137776 1.246423 1.035300 1.277320  1.187472
# 100% 12.329442 17.515992 7.654660     4.900264 8.000668 9.451632 9.283560 9.712032 9.719232 8.319498 31.475180
#           drop  Excit_01 Inhib_04 Excit_08 Excit_14 Excit_02 Excit_06 Excit_07
# 0%    0.000000  0.023800 0.143892 0.019960 0.081284 0.016268 0.199600 0.137636
# 25%   0.213146  0.476444 0.345158 0.321237 0.477489 0.469060 0.356304 1.230960
# 50%   0.732060  0.742440 0.759984 0.778440 0.864066 0.902736 1.057420 1.477152
# 75%   1.773100  1.166396 1.447852 1.533816 1.438500 1.387359 1.712640 1.785096
# 100% 12.068628 13.872384 9.825768 8.841504 6.457070 8.640852 5.808176 4.430430

sapply(cellType.idx, function(x) {
  quantile(sce$sum[x])
})[, order(sapply(cellType.idx, function(x) {
  quantile(sce$sum[x])["50%"]
}))]
#      Excit_15  drop Excit_09 Oligo_03 Oligo_01 Micro EndoMural_02   Astro Excit_13 Oligo_02 Excit_10 Excit_05      OPC
# 0%        678   220      462      359   300.00   350        407.0   276.0    648.0    968.0      533   292.00   397.00
# 25%      1101  1214     1389     1464  1592.25  1756       2020.5  1474.0   2447.5   3750.0     3298  3077.25  3933.25
# 50%      1218  1856     2121     2211  2554.00  2633       3140.0  3280.0   3445.0   5113.5     5489  5668.50  5951.00
# 75%      1348  3067     3225     3410  3879.25  4300       4546.0  7305.5   4680.5   7314.5     8100 10000.00  8510.00
# 100%    15127 27828    50703    36108 16212.00 15576      32302.0 61842.0  16113.0  49758.0    95913 65772.00 83867.00
#      Excit_14 EndoMural_01 Inhib_03 Inhib_06  Inhib_05 Inhib_01 Inhib_04 Excit_11 Excit_12 Excit_03 Inhib_02  Excit_07
# 0%    2853.00      2040.00   641.00   1598.0   1402.00   317.00     1571  4585.00  1215.00     1651   4438.0   2691.00
# 25%   5537.75      6102.75  7002.00   8979.5  10839.25 10929.00    12837 14424.25 13997.75    14598  16177.5  17929.75
# 50%   7893.00      8219.50 12173.00  13173.0  15142.50 16488.00    18805 19989.50 20810.00    20879  23964.0  25357.50
# 75%  17635.75     11667.25 18893.25  20914.5  21695.25 23337.75    26787 28125.75 28306.25    30108  35836.5  39349.00
# 100% 74260.00     44112.00 61063.00  99561.0 135464.00 97656.00    75908 67166.00 59560.00   103322 118342.0 105390.00
#      Excit_04 Excit_01 Excit_08 Excit_06 Excit_02
# 0%     1449.0    327.0    922.0     1437    629.0
# 25%   21769.5  22253.5  19554.5    24773  27708.5
# 50%   31522.0  33633.0  33928.0    37075  42935.0
# 75%   47549.5  51601.5  53335.0    51649  65244.0
# 100% 179194.0 203794.0 188090.0   152333 296099.0

sapply(cellType.idx, function(x) {
  quantile(sce$subsets_Mito_percent[x])
})[, order(sapply(cellType.idx, function(x) {
  quantile(sce$subsets_Mito_percent[x])["50%"]
}))]
#         Excit_06    Excit_08    Excit_02   Excit_07   Excit_01   Inhib_04   Excit_03   Inhib_05    Excit_04  Inhib_06
# 0%   0.002438608 0.001886792 0.002068851 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.003583716 0.0000000
# 25%  0.042550605 0.056380879 0.065285798 0.06414893 0.06867666 0.08194802 0.07818608 0.08390369 0.117383520 0.1139963
# 50%  0.084033613 0.109454614 0.122737036 0.12342525 0.12580109 0.13477997 0.13777358 0.14170341 0.175478065 0.1833181
# 75%  0.180879403 0.222722285 0.244560797 0.22520799 0.22560217 0.25029052 0.25425884 0.23491904 0.312303374 0.2918794
# 100% 3.612059158 5.267234702 4.861693210 2.23619516 4.47400242 3.28638498 3.87643852 3.16921336 2.887537994 3.5264901
#       Inhib_02       OPC   Excit_12   Excit_11  Inhib_03  Oligo_02  Inhib_01 EndoMural_01     Astro     Micro
# 0%   0.0000000 0.0000000 0.01966491 0.01526718 0.0000000 0.0000000 0.0000000    0.1190090 0.0000000 0.0000000
# 25%  0.1010100 0.1244257 0.13335564 0.14361395 0.1620331 0.1866169 0.2183281    0.3309703 0.2677177 0.3186743
# 50%  0.1917883 0.2128511 0.23253048 0.24706812 0.2986780 0.3133674 0.4091633    0.4636469 0.5033557 0.5308392
# 75%  0.3264469 0.3629894 0.34917901 0.37872377 0.6385777 0.5284328 0.7746272    0.6182482 1.0065139 0.8727477
# 100% 3.0198447 4.4997040 2.85996055 3.66186919 7.0487994 5.5207315 4.7677262    2.6216485 7.6875387 5.6744186
#      EndoMural_02  Excit_14  Oligo_03   Excit_10  Excit_05  Excit_13  Excit_09  Excit_15  Oligo_01      drop
# 0%      0.0000000 0.0168672 0.0000000 0.06800408 0.0000000 0.2733663 0.0936427 0.2906977  0.000000  0.000000
# 25%     0.5050875 0.3041596 0.3949643 0.62885069 0.4435092 1.4797343 1.0055866 1.4754156  2.003345  2.858362
# 50%     0.7608696 0.8055773 0.8540014 0.88545592 0.9957153 1.8018018 1.8310691 1.8327272  3.086420  4.385307
# 75%     1.4159816 1.0237141 2.6376888 1.50116521 2.1279506 2.2319259 2.9976019 2.2038383  4.387429  5.580677
# 100%    7.0837167 2.7013507 8.9666951 7.51043115 8.0777538 6.1119293 8.2758621 3.9635355 10.294118 12.380192


#### Check out sub-clusters of Excit_13 with Astro-like expression
sce_Excit13 <- sce[, sce$cellType_hc == "Excit_13"]

sce_Excit13$prelimCluster <- droplevels(sce_Excit13$prelimCluster)

table(sce_Excit13$BrNum, sce_Excit13$prelimCluster)
#         28 133 149 210 236 255 291 295
# Br2720  12   7  98   3   0   0   0  41
# Br2743   2   5   5   0   0   0   0   0
# Br3942 347   2  43   0   0   1   0   0
# Br6423  11 422  72  53 157   1   0   0
# Br6432   0   0   1   0   0   0   0   0
# Br6471   1   2   1   0   0   0   0   0
# Br6522   2   2   3   0   0   0   0   0
# Br8492  73   2  25   0   0   0   0   0
# Br8667  14   4   1   2   5  26 121   0

table(sce_Excit13$Sample[sce_Excit13$prelimCluster %in% c(133, 210, 236)])
table(sce_Excit13$subject[sce_Excit13$prelimCluster %in% c(133, 210, 236)])

table(unfactor(sce$prelimCluster)[sce$cellType_hc == "Excit_13"])
# 133 149 210 236 255  28 291 295
# 446 249  58 162  28 462 121  41

my_plotMarkers(
  sce = sce[, sce$cellType_hc == "Excit_13"],
  marker_list = markers.mathys.tran[c("astrocyte", "excit_neuron")],
  cat = "prelimCluster",
  pdf_fn = here(plot_dir, "markers_Excit_13_check.pdf")
)

dbs_excit13 <- ggcells(sce[, sce$cellType_hc == "Excit_13"],
                       mapping = aes(x = prelimCluster, y = doubletScore, fill = prelimCluster)
) +
  geom_boxplot() +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "None", axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

ggsave(dbs_excit13, filename = here(plot_dir, "qc_Excit_13_doublet_scores.png"))

# #### Explore Prelim clusters ####
# prelim_mito <- ggcells(sce, mapping = aes(
#     x = reorder(prelimCluster, subsets_Mito_percent, FUN = median),
#     y = subsets_Mito_percent, color = cellType_hc
# )) +
#     geom_boxplot() +
#     scale_color_manual(values = cell_type_colors) +
#     theme_bw() +
#     theme(
#         legend.position = "None", axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1)
#     )
#
# ggsave(prelim_mito, filename = here(plot_dir, "prelim_mitoRate.png"), width = 30)
#
#
# prelim_UMI <- ggcells(sce, mapping = aes(
#     x = reorder(prelimCluster, sum, FUN = median),
#     y = sum, color = cellType_hc
# )) +
#     geom_boxplot() +
#     scale_color_manual(values = cell_type_colors) +
#     theme_bw() +
#     theme(
#         legend.position = "None", axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1)
#     )
#
# ggsave(mito_ct_hc, filename = here(plot_dir, "prelim_UMI.png"), width = 30)
#
# prelim_doublets <- ggcells(sce, mapping = aes(
#     x = reorder(prelimCluster, doubletScore, FUN = median),
#     y = doubletScore, color = cellType_hc
# )) +
#     geom_boxplot() +
#     scale_color_manual(values = cell_type_colors) +
#     theme_bw() +
#     theme(
#         legend.position = "None", axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1)
#     )
#
# ggsave(prelim_doublets, filename = here(plot_dir, "prelim_doublets.png"), width = 30)

sgejobs::job_single('09_cluster_annotation_explore', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 09_cluster_annotation_explore.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
