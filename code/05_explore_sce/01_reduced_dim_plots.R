library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("here")

#### Load Data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

#### Set Up Plotting ####
my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "05_explore_sce", "01_reduced_dim_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]

#### UMAP plots ####
UMAP_cellTypes_hc <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors) +
    my_theme +
    coord_equal()

ggsave(UMAP_cellTypes_hc +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_cellType.png"), width = 9
)

## pdf versions
ggsave(UMAP_cellTypes_hc +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_cellType.pdf"), width = 9
)

ggsave(UMAP_cellTypes_hc +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_no_legend.png")
)

## Add facet for cell types
ggsave(UMAP_cellTypes_hc +
    facet_wrap(~cellType_hc) +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_facet.png"), width = 10, height = 10
)

## layer annotaions
UMAP_cellType_layer <- ggcells(sce[, !is.na(sce$cellType_layer)], mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType_layer)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors_layer) +
    my_theme +
    coord_equal()

ggsave(UMAP_cellType_layer +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_cellType_layer.png"), width = 9
)

ggsave(UMAP_cellType_layer +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_layer_no_legend.png")
)


#### Plot clusters in TSNE ####
TSNE_HC_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors) +
    my_theme +
    coord_equal()

ggsave(TSNE_HC_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "TSNE_cellType.png"), width = 10
)

## Add facet for cell types
ggsave(TSNE_HC_cellTypes +
    facet_wrap(~cellType_hc) +
    theme(legend.position = "None"),
filename = here(plot_dir, "TSNE_cellType_facet.png"), width = 10, height = 10
)


TSNE_broad_cellTypes <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2)) +
    geom_point(aes(colour = cellType_broad_hc), size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors[levels(sce$cellType_broad_hc)], drop = TRUE) +
    my_theme +
    coord_equal()

ggsave(TSNE_broad_cellTypes +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "TSNE_broad_cellType.png"), width = 10
)

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


UMI_ct_hc <- ggcells(sce, mapping = aes(x = cellType_hc, y = sum, fill = cellType_hc)) +
    geom_boxplot() +
    scale_fill_manual(values = cell_type_colors) +
    my_theme +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave(UMI_ct_hc, filename = here(plot_dir, "UMI_HC_cellType.png"), width = 10)


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


dbs_ct_hc <- ggcells(sce, mapping = aes(x = cellType_hc, y = doubletScore, fill = cellType_hc)) +
    geom_boxplot() +
    scale_fill_manual(values = cell_type_colors) +
    my_theme +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

ggsave(dbs_ct_hc, filename = here(plot_dir, "doubletScore_HC_cellType.png"), width = 10)


cellType.idx <- splitit(sce$cellType_hc)
# sapply(c("Excit", "Inhib", "MSN"), function(x){grep(x, names(cellType.idx))})

sapply(cellType.idx, function(x) {
    quantile(sce$doubletScore[x])
})[, order(sapply(cellType.idx, function(x) {
    quantile(sce$doubletScore[x])["50%"]
}))]
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

sapply(cellType.idx, function(x) {
    quantile(sce$sum[x])
})[, order(sapply(cellType.idx, function(x) {
    quantile(sce$sum[x])["50%"]
}))]
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

dat_symb <- rownames(sce)[match(rownames(dat_small), rowData(sce)$gene_id)]
any(is.na(dat_symb))

layer_idx <- splitit(dat_small$Layer)
layer_top4 <- purrr::map(layer_idx, ~ dat_symb[.x][1:4])
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
my_plotMarkers(
    sce = sce,
    marker_list = layer_top4,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "HC_layer_markers_ct.pdf")
)

my_plotMarkers(
    sce = sce,
    marker_list = markers.mathys.tran,
    cat = "cellType_k",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "mb_kmeans_29_mathys_markers_ct.pdf")
)

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
