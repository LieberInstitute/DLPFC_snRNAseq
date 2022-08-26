library("SingleCellExperiment")
library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")
library("numform")
library("jaffelab")

source("utils.R")
source("my_plotExpression.R")
## load data
## TODO USE HDF5
load(here("processed-data", "03_build_sce", "sce_harmony_Sample.Rdata"), verbose = TRUE)
load(here("processed-data", "03_build_sce", "km_res.Rdata"), verbose = TRUE)

#### Define SCE clusters with chosen k ####
## picked k = 29
k_list <- seq(5, 50) ## keep index
sce$kmeans <- as.factor(paste0("mbk", f_pad_zero(km_res[[which(k_list == 29)]]$Clusters)))

(cluster_tab <- table(sce$kmeans))
# mbk01 mbk02 mbk03 mbk04 mbk05 mbk06 mbk07 mbk08 mbk09 mbk10 mbk11 mbk12 mbk13 mbk14 mbk15 mbk16 mbk17 mbk18 mbk19 mbk20
# 4233   264  1420 28963    35     2  2791   733  1235  1146  1192    16  4659   174  1587  3542   466  1645  1518   127
# mbk21 mbk22 mbk23 mbk24 mbk25 mbk26 mbk27 mbk28 mbk29
# 536  2718  1427  1615  5900   113  1193  2235  6119

summary(as.numeric(table(sce$kmeans)))

table(sce$kmeans, sce$round)
table(sce$kmeans, sce$Sample)
table(sce$kmeans, sce$subject)

#### check doublet score for each prelim clust ####
clusIndexes <- splitit(sce$kmeans)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii) {
    median(sce$doubletScore[ii])
})

## pretty low doublet scores
round(prelimCluster.medianDoublet, 2)
# mbk01 mbk02 mbk03 mbk04 mbk05 mbk06 mbk07 mbk08 mbk09 mbk10 mbk11 mbk12 mbk13 mbk14 mbk15 mbk16 mbk17 mbk18 mbk19 mbk20
# 0.41  0.40  0.70  0.54  1.26  0.55  0.16  0.77  0.72  1.01  0.20  0.78  0.46  0.82  0.56  0.14  0.60  0.86  0.80  0.77
# mbk21 mbk22 mbk23 mbk24 mbk25 mbk26 mbk27 mbk28 mbk29
# 0.15  0.38  0.67  1.58  0.67  0.40  1.49  0.74  1.13
summary(prelimCluster.medianDoublet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1370  0.4148  0.6658  0.6797  0.8024  1.5806

#### Plot TSNE ####
# Plotting set up
load(here("processed-data", "03_build_sce", "color_palletes.Rdata"), verbose = TRUE)
names(iWantHue_k29) <- levels(sce$kmeans)
cluster_colors <- iWantHue_k29
# cluster_colors <- DeconvoBuddies::create_cell_colors(cell_types = levels(sce$kmeans), pallet = "gg")

my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "03_build_sce", "cluster")

## Plot clusters in TSNE
TSNE_clusters <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = kmeans)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cluster_colors) +
    my_theme +
    coord_equal()

ggsave(TSNE_clusters +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 2))),
filename = here(plot_dir, "clusters_mbkm-29.png"), width = 10
)

ggsave(TSNE_clusters +
    facet_wrap(~kmeans) +
    theme(legend.position = "none"),
filename = here(plot_dir, "clusters_mbkm-29_facet.png"), width = 10, height = 10
)

#### Marker Genes ####
# Just for logcounts
sce <- batchelor::multiBatchNorm(sce, batch = sce$Sample)

# load Mathy's markers
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/markers.rda", verbose = TRUE)
# markers.mathys.custom
all(unlist(markers.mathys.custom) %in% rowData(sce)$gene_name)

rownames(sce) <- rowData(sce)$gene_name

source("my_plotExpression.R")

pdf(here("plots", "03_build_sce", "cluster", "mb_kmeans_29_mathys_markers.pdf"), height = 6, width = 8)
for (i in 1:length(markers.mathys.custom)) {
    message(names(markers.mathys.custom)[[i]])
    print(
        my_plotExpression(sce,
            genes = markers.mathys.custom[[i]],
            title = names(markers.mathys.custom)[[i]],
            cat = "kmeans",
            fill_colors = cluster_colors
        )
    )
}
dev.off()

#### Tran Maynard Top Markers####
dlpfc_markers <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv")
dlpfc_markers_list <- as.list(dlpfc_markers[1:4, grepl("_1vAll", colnames(dlpfc_markers))])

pdf(here("plots", "03_build_sce", "cluster", "mb_kmeans_29_Tran_markers.pdf"), height = 6, width = 8)
for (i in 1:length(dlpfc_markers_list)) {
    message(names(dlpfc_markers_list)[[i]])
    f <- dlpfc_markers_list[[i]]
    f_good <- f[f %in% rownames(sce)]
    if (f != f_good) message("Missing...", paste(f[!f %in% f_good], collapse = ", "))
    print(
        my_plotExpression(sce,
            genes = f_good,
            title = names(dlpfc_markers_list)[[i]],
            cat = "kmeans",
            fill_colors = cluster_colors
        )
    )
}
dev.off()

#### Mean Ratio Top Markers####
load("/dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

mr_list <- marker_stats %>%
    ungroup() %>%
    arrange(cellType.target) %>%
    filter(rank_ratio < 5) %>%
    select(rank_ratio, cellType.target, Symbol) %>%
    tidyr::pivot_wider(names_from = "cellType.target", values_from = "Symbol") %>%
    select(-rank_ratio) %>%
    as.list()

pdf(here("plots", "03_build_sce", "cluster", "mb_kmeans_29_MeanRatio_markers.pdf"), height = 6, width = 8)
for (i in 1:length(mr_list)) {
    message(names(mr_list)[[i]])
    f <- mr_list[[i]]
    f_good <- f[f %in% rownames(sce)]
    if (f != f_good) message("Missing...", paste(f[!f %in% f_good], collapse = ", "))
    print(
        my_plotExpression(sce,
            genes = f_good,
            title = names(mr_list)[[i]],
            cat = "kmeans",
            fill_colors = cluster_colors
        )
    )
}
dev.off()

#### Annotate with marker cell types ####
anno <- read.csv(here("processed-data", "03_build_sce", "DLPFC_k29_anno.csv"))
table(anno$broad)

sce$cellType.broad <- factor(anno$broad[match(sce$kmeans, anno$cluster)])
table(sce$cellType.broad)
# Astro Excit Inhib Micro Nural Oligo   OPC small
# 3557 18035 11380  2242  1330 39257  1791    12

## Precentage cell composition
(prop <- 100 * round(table(sce$cellType.broad) / ncol(sce), 3))
# Astro Excit Inhib Micro Nural Oligo   OPC small
# 4.6  23.2  14.7   2.9   1.7  50.6   2.3   0.0

## plot broad
load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/processed-data/00_data_prep/cell_colors.Rdata", verbose = TRUE)

pdf(here("plots", "03_build_sce", "cluster", "mb_broad_mathys_markers.pdf"), height = 6, width = 8)
for (i in 1:length(markers.mathys.custom)) {
    message(names(markers.mathys.custom)[[i]])
    print(
        plotExpressionCustom(
            sce = sce,
            features = markers.mathys.custom[[i]],
            features_name = names(markers.mathys.custom)[[i]],
            anno_name = "cellType.broad"
        ) +
            scale_color_manual(values = cell_colors)
    )
}
dev.off()

save(sce, here("processed-data", "sce", "sce_DLPFC.Rdata"))

#### TRan Maynard data check ####
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda", verbose = TRUE)
t <- table(ss(as.character(sce.dlpfc$cellType), "_"))
# Astro      Excit      Inhib Macrophage      Micro      Mural      Oligo        OPC      Tcell
# 782       2388       1580         10        388         18       5455        572          9
(prop_tran <- 100 * round(t / ncol(sce.dlpfc), 3))
# Astro      Excit      Inhib Macrophage      Micro      Mural      Oligo        OPC      Tcell
# 7.0       21.3       14.1        0.1        3.5        0.2       48.7        5.1        0.1

data.frame(prop)

message("Saving Data - ", Sys.time())
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce", "sce_DLPFC"))
