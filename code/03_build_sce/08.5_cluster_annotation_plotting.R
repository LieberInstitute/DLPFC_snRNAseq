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

plot_dir <- here("plots", "03_build_sce", "08.5_cluster_annotation_plotting")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Load Data ####
load(here("processed-data", "sce","sce_DLPFC.Rdata"), verbose = TRUE)

#### Replot with Annotations ####
all(ct %in% names(cell_type_colors)) ## if not revisit cell_colors

## Plot clusters in TSNE
TSNE_HC_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(TSNE_HC_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "TSNE_HC-29_cellType.png"), width = 10
)


TSNE_HC_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2)) +
    geom_point(aes(colour = cellType_broad_hc), size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_broad_hc)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(TSNE_HC_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "TSNE_HC_broad_cellType.png"), width = 10
)


## UMAP
UMAP_HC_cellTypes <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(UMAP_HC_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_HC_cellType.png"), width = 10
)


UMAP_HC_cellTypes <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2)) +
    geom_point(aes(colour = cellType_broad_hc), size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_broad_hc)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(UMAP_HC_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_HC_broad_cellType.png"), width = 10
)

TSNE_km_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_k)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_k)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(TSNE_km_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "TSNE_mbkm-29_cellType.png"), width = 10
)


## markers
my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "HC_mathys_markers_ct.pdf")
)

my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "cellType_k",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "mb_kmeans_29_mathys_markers_ct.pdf")
)


#### Compare annotations ####
table(sce$kmeans)
table(sce$collapsedCluster)

cluster_compare <- table(sce$kmeans, sce$collapsedCluster)
cluster_compare_prop <- sweep(cluster_compare, 2, colSums(cluster_compare), `/`)

cluster_compare > 0

hc_counts <- table(sce$kmeans)

hc_anno <- anno_hc %>%
    select(cluster, broad) %>%
    tibble::column_to_rownames("cluster") %>%
    mutate(n = table(sce$collapsedCluster))
km_anno <- anno_k %>%
    select(cluster, broad) %>%
    tibble::column_to_rownames("cluster") %>%
    mutate(n = table(sce$kmeans))

broad_pallet <- cell_type_colors_broad[unique(c(hc_anno$broad, km_anno$broad))]

png(here(plot_dir, "cluster_compare_heatmap.png"), height = 800, width = 800)
pheatmap(cluster_compare_prop,
    annotation_col = hc_anno,
    annotation_row = km_anno,
    annotation_colors = list(broad = temp_pallet)
)
dev.off()

## Jaccard indicies
library(bluster)
jacc.mat <- linkClustersMatrix(sce$kmeans, sce$collapsedCluster)

png(here(plot_dir, "cluster_compare_heatmap_jacc.png"), height = 800, width = 800)
pheatmap(jacc.mat,
    annotation_col = hc_anno,
    annotation_row = km_anno,
    annotation_colors = list(broad = temp_pallet)
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

hc_ct_anno <- table(sce$cellType_hc) %>%
    as.data.frame() %>%
    rename(ct = Var1, n = Freq)
rownames(hc_ct_anno) <- hc_ct_anno$ct

km_ct_anno <- table(sce$cellType_k) %>%
    as.data.frame() %>%
    rename(ct = Var1, n = Freq)
rownames(km_ct_anno) <- km_ct_anno$ct

png(here(plot_dir, "cellType_compare_heatmap.png"), height = 800, width = 800)
pheatmap(ct_tab_prop,
    annotation_col = hc_ct_anno,
    annotation_row = km_ct_anno,
    annotation_colors = list(ct = cell_type_colors)
)
dev.off()

jacc.mat.ct <- linkClustersMatrix(sce$cellType_k, sce$cellType_hc)

png(here(plot_dir, "cellType_compare_heatmap_jacc.png"), height = 800, width = 800)
pheatmap(jacc.mat.ct,
    annotation_col = hc_ct_anno,
    annotation_row = km_ct_anno,
    annotation_colors = list(ct = cell_type_colors)
)
dev.off()

## compare broad cell types

(ctb_tab <- table(sce$cellType_broad_k, sce$cellType_broad_hc))
ctb_tab_prop <- sweep(ctb_tab, 2, colSums(ctb_tab), `/`)

hc_broad_anno <- table(sce$cellType_broad_hc) %>%
    as.data.frame() %>%
    rename(broad = Var1, n = Freq)
rownames(hc_broad_anno) <- hc_broad_anno$broad

km_broad_anno <- table(sce$cellType_broad_k) %>%
    as.data.frame() %>%
    rename(broad = Var1, n = Freq)
rownames(km_broad_anno) <- km_broad_anno$broad

png(here(plot_dir, "cellType_broad_compare_heatmap.png"), height = 800, width = 800)
pheatmap(ctb_tab_prop,
    annotation_col = hc_broad_anno,
    annotation_row = km_broad_anno,
    annotation_colors = list(broad = broad_pallet)
)
dev.off()

jacc.mat.broad <- linkClustersMatrix(sce$cellType_broad_k, sce$cellType_broad_hc)

png(here(plot_dir, "cellType_broad_compare_heatmap_jacc.png"), height = 800, width = 800)
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

ggsave(UMI_ct_k, filename = here(plot_dir, "UMI_mbkm-29_cellType.png"), width = 10)


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

ggsave(dbs_ct_k, filename = here(plot_dir, "doubletScore_mbkm-29_cellType.png"), width = 10)


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


cellType.idx <- splitit(sce$cellType_hc)
# sapply(c("Excit", "Inhib", "MSN"), function(x){grep(x, names(cellType.idx))})

sapply(cellType.idx, function(x) {
    quantile(sce$doubletScore[x])
})[, order(sapply(cellType.idx, function(x) {
    quantile(sce$doubletScore[x])["50%"]
}))]
#       Excit_15    Micro  Oligo_03 Excit_09     Astro       OPC Excit_10   Oligo_02 Excit_13  Excit_11 EndoMural_02
# 0%   0.0086180 0.000000  0.000000  0.00000  0.000000  0.000000 0.000000  0.0000000 0.000000  0.050708      0.000000
# 25%  0.0348360 0.032112  0.042816  0.04134  0.057876  0.068640 0.077310  0.0811380 0.103416  0.144880      0.085632
# 50%  0.0751140 0.099216  0.137040  0.13744  0.143964  0.166296 0.180978  0.1896200 0.197182  0.282516      0.300958
# 75%  0.1874415 0.388596  0.568988  0.45220  0.345052  0.639158 0.426949  0.4642535 0.426949  0.476293      0.765917
# 100% 4.4314060 6.411696 15.852624  6.83060 12.939216 10.329944 7.175700 16.9815360 8.621844 24.470232      8.387922
#       Inhib_03  Excit_12 Inhib_02 EndoMural_01 Excit_03 Inhib_06 Inhib_05 Excit_05 Excit_04  Inhib_01 Oligo_01
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

sapply(cellType.idx, function(x) {
    quantile(sce$sum[x])
})[, order(sapply(cellType.idx, function(x) {
    quantile(sce$sum[x])["50%"]
}))]
#      Excit_15 Oligo_01 Excit_09 Oligo_03 Micro EndoMural_02   Astro Excit_13 Oligo_02 Excit_10 Excit_05      OPC
# 0%        678      220      462      359   350         407.0   276.0    648.0    968.0      533   292.00   397.00
# 25%      1101     1232     1389     1464  1756        2020.5  1474.0   2447.5   3750.0     3298  3077.25  3933.25
# 50%      1218     1907     2121     2211  2633        3140.0  3280.0   3445.0   5113.5     5489  5668.50  5951.00
# 75%      1348     3147     3225     3410  4300        4546.0  7305.5   4680.5   7314.5     8100 10000.00  8510.00
# 100%    15127    27828    50703    36108 15576       32302.0 61842.0  16113.0  49758.0    95913 65772.00 83867.00
#      Excit_14 EndoMural_01 Inhib_03 Inhib_06  Inhib_05 Inhib_01 Inhib_04 Excit_11 Excit_12 Excit_03 Inhib_02  Excit_07
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


#### ALL QC plots ####
library(patchwork)

ggsave((UMI_ct_hc / dbs_ct_hc / mito_ct_hc), filename = here(plot_dir, "qc_ALL_HC_cellType.png"), height = 7.5)


## Mito over UMAP
library(scales)
UMAP_mito <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2)) +
    geom_point(aes(colour = subsets_Mito_percent), size = 0.2, alpha = 0.3) +
    scale_color_continuous(type = "viridis") +
    my_theme +
    coord_equal()

ggsave(UMAP_mito, filename = here(plot_dir, "UMAP_Mito.png"), width = 10)



#### Explore layer markers ####
## may look at cluster marker genes instead
# load("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/07_spatial_registration/t_cor_plot_top_genes_k7.rda", verbose = TRUE)
# # dat_small
# table(dat_small$Layer)
# # Layer1 Layer2 Layer3 Layer4 Layer5 Layer6     WM
# # 91     87     65     80     78     83    100
#
# dat_symb <- rownames(sce)[match(rownames(dat_small),rowData(sce)$gene_id)]
# any(is.na(dat_symb))
#
# layer_idx <- splitit(dat_small$Layer)
# layer_top4 <- purrr::map(layer_idx, ~dat_symb[.x][1:4])


## markers
my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "cellType_k",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_mathys_ct_mb_kmeans_29.pdf")
)

my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_mathys_ct_hc_29.pdf")
)

# my_plotMarkers(sce = sce,
#                marker_list = layer_top4,
#                cat = "cellType_hc",
#                fill_colors = cell_type_colors,
#                pdf_fn = here(plot_dir, "HC_layer_markers_ct.pdf"))


## Inhib sub-types
markers.Inhib_subtypes <- list(
    "inhib_markers pg1" = c("HTR3A", "VIP", "CCK", "NPY", "CRHBP", "CALB2", "PNOC"),
    "inhib_Zhaung2022" = c(
        "SP8", # VIP
        "KLF5", # SST
        "LGI2", # PVALB
        "LAMP5"
    ) # LAMP5
)

my_plotMarkers(
    sce = sce[, grep("Inhib", sce$cellType_hc)],
    marker_list = markers.Inhib_subtypes,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "HC_inhib_markers_ct.pdf")
)

sce_inhib <- sce[, grep("Inhib", sce$cellType_hc)]
sce_inhib$cellType_hc <- droplevels(sce_inhib$cellType_hc)
sce_inhib$collapsedCluster <- droplevels(sce_inhib$collapsedCluster)
table(sce_inhib$cellType_hc, sce_inhib$collapsedCluster)


#### Check out sub-clusters of Excit_13 with Astro-like expression
sce_Excit13 <- sce[, sce$cellType_hc == "Excit_13"]

sce_Excit13$prelimCluster <- droplevels(sce_Excit13$prelimCluster)

table(sce_Excit13$subject, sce_Excit13$prelimCluster)
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

prelimCluster.medianDoublet[c(28, 133, 149, 210, 236, 255, 291, 295)]

## questionable
prelimCluster.medianDoublet[c(133, 210, 236)]

my_plotMarkers(
    sce = sce[, sce$cellType_hc == "Excit_13"],
    marker_list = markers.mathys.tran[c("astrocyte", "excit_neuron")],
    cat = "prelimCluster",
    pdf_fn = here(plot_dir, "Excit_13_check.pdf")
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

ggsave(dbs_excit13, filename = here(plot_dir, "Excit_13_doublet_scores.png"))


#### Explore Prelim clusters ####

prelim_mito <- ggcells(sce, mapping = aes(
    x = reorder(prelimCluster, subsets_Mito_percent, FUN = median),
    y = subsets_Mito_percent, color = cellType_hc
)) +
    geom_boxplot() +
    scale_color_manual(values = cell_type_colors) +
    theme_bw() +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave(prelim_mito, filename = here(plot_dir, "prelim_mitoRate.png"), width = 30)


prelim_UMI <- ggcells(sce, mapping = aes(
    x = reorder(prelimCluster, sum, FUN = median),
    y = sum, color = cellType_hc
)) +
    geom_boxplot() +
    scale_color_manual(values = cell_type_colors) +
    theme_bw() +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave(mito_ct_hc, filename = here(plot_dir, "prelim_UMI.png"), width = 30)

prelim_doublets <- ggcells(sce, mapping = aes(
    x = reorder(prelimCluster, doubletScore, FUN = median),
    y = doubletScore, color = cellType_hc
)) +
    geom_boxplot() +
    scale_color_manual(values = cell_type_colors) +
    theme_bw() +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave(prelim_doublets, filename = here(plot_dir, "prelim_doublets.png"), width = 30)
