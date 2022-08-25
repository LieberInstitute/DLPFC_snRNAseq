library("SingleCellExperiment")
# library("SeuratData")
library("patchwork")
library("here")
library("pheatmap")
library("scater")
library("bluster")
library("tidyverse")
library("sessioninfo")
library(dplyr)

#### Plot Setup ####
plot_dir = here("plots","05_explore_sce","06_explore_azimuth_annotations")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"), verbose = TRUE)


source(here("code", "03_build_sce","my_plotExpression.R"))

## load data
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

pd <- as.data.frame(colData(sce))

cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]
cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]

#### prop breakdown ####

ct_counts <- pd |> 
  group_by(cellType_hc, cellType_azimuth) |>
  summarize(n_cell = n())

azimuth_prop <- pd |>
  count(cellType_azimuth) |> 
  left_join(ct_counts) |>
  mutate(prop = n_cell/n)

azimuth_prop_bar <- azimuth_prop |>
  ggplot(aes(x = cellType_azimuth, y = prop, fill = cellType_hc)) +
  geom_col() +
  scale_fill_manual(values = cell_type_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(azimuth_prop_bar, filename = here(plot_dir, "prop_bar_azimuth.png"))

azimuth_count_bar <- azimuth_prop |>
  ggplot(aes(x = cellType_azimuth, y = n_cell, fill = cellType_hc)) +
  geom_col() +
  scale_fill_manual(values = cell_type_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(azimuth_count_bar, filename = here(plot_dir, "count_bar_azimuth.png"))

library(patchwork)
ggsave(azimuth_prop_bar + theme(legend.position = "None") + azimuth_count_bar, filename = here(plot_dir, "count-prop_bar_azimuth.png"), width = 12)

#### Azimuth Cell Types ####

az <- unique(sce$cellType_azimuth)
gsub("(L[1-6]).*", "\\1", az)
# [1] "Sst"       "Oligo"     "OPC"       "Vip"       "L6"        "L5"        "VLMC"      "Astro"     "Endo"     
# [10] "L2"        "Micro-PVM" "L6"        "Sst Chodl" "Lamp5"     "Pvalb"     "Sncg"      "L6"        "L6"       
# [19] "L5"        "L5"

azimuth_cellType_notes <- data.frame(azimuth = az,
                                     cellType_broad = c("Inhib","Oligo", "OPC", "Inhib","Excit", 
                                                   "Excit", "EndoMural","Astro","EndoMural",
                                                   "Excit", "Micro", "Excit", "Excit", "Inhib", 
                                                   "Inhib","Excit","Excit","Excit","Excit","Excit")) |>
  mutate(Layer = ifelse(grepl("L[1-6]", azimuth), gsub("(L[1-6]|L[1-6]/[1-6]).*", "\\1", azimuth), NA))


#### Compare to HC annotations 
ct_annos <- table(sce$cellType_hc, sce$cellType_azimuth)

png(here(plot_dir, "azimuth_v_hc.png"),height = 800, width = 800)
pheatmap(ct_annos,
         main = "Number of Nuclei",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

## Prop of HC annotation groups
ct_annos_prop <- sweep(ct_annos,1,rowSums(ct_annos),`/`)

ct_max <- apply(ct_annos_prop, 1, max)

png(here(plot_dir, "azimuth_v_hc-prop.png"),height = 800, width = 800)
pheatmap(ct_annos_prop,
         main = "Proportion of HC Cluster",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

png(here(plot_dir, "azimuth_v_hc-prop-cluster.png"),height = 800, width = 800)
pheatmap(ct_annos_prop,
         main = "Number of Nuclei")
dev.off()

## How do clusters correspond?
jacc.mat <- linkClustersMatrix(sce$cellType_hc, sce$cellType_azimuth)

png(here(plot_dir, "azimuth_v_hc-jacc.png"),height = 800, width = 800)
pheatmap(jacc.mat,
         main = "Strength of the correspondence between Azimuth & HC")
dev.off()


#### Annotate the Heatmap ####
library(viridis)
jacc.mat[1:5, 1:5]

## Number of markers
load(here("processed-data", "03_build_sce","cell_type_markers.Rdata"), verbose = TRUE)
n_markers <- sapply(markers_1vALL, function(x) sum(x$FDR < 0.05))

## Spatial Registration
load(here("processed-data","05_explore_sce","07_sptatial_registration","layer_cor.Rdata"), verbose = TRUE)
max_cor <- apply(cor, 1, max)
max_layer <- colnames(cor)[apply(cor,1,which.max)[names(n_markers)]]

## build annotation df
hc_anno <- data.frame(cellType_broad = gsub("\\.|_[0-9]+","",names(n_markers)),
                      n_markers = n_markers,
                      max_cor = max_cor[names(n_markers)],
                      Layer = gsub("ayer","",max_layer))

rownames(hc_anno) <- gsub("\\.","",rownames(hc_anno))


azimuth_anno <- azimuth_cellType_notes |> 
  column_to_rownames("azimuth")

layers <- c("WM", "L1", "L2", "L2/3", "L3","L4","L5","L5/6","L6")

azimuth_anno$Layer <- factor(azimuth_anno$Layer, levels = layers)
hc_anno$Layer <- factor(hc_anno$Layer, levels = layers)

layer_colors <- viridis(length(layers))
names(layer_colors) <- layers 

png(here(plot_dir, "azimuth_v_hc_annotation.png"),height = 800, width = 800)
pheatmap(jacc.mat,
         color= inferno(100),
         annotation_row = hc_anno,
         annotation_col = azimuth_anno,
         annotation_colors = list(Layer = layer_colors, cellType_broad = cell_type_colors_broad),
         main = "Strength of the correspondence between Azimuth & HC")
dev.off()

## Prop of prelim annotations

prelim_anno <- pd |> 
  group_by(prelimCluster, cellType_hc) |> 
  count() |>
  column_to_rownames("prelimCluster") |>
  arrange(cellType_hc)

prelim_anno_n <- table(sce$prelimCluster, sce$cellType_hc)
prelim_anno_n <- table(sce$prelimCluster, sce$cellType_azimuth)
prelim_anno_prop <- sweep(prelim_anno_n,1,rowSums(prelim_anno_n),`/`)

table(sce$prelimCluster, sce$cellType_hc)
ct_max <- apply(ct_annos_prop, 1, max)


png(here(plot_dir, "azimuth_v_prelim.png"),height = 1200, width = 800)
pheatmap(prelim_anno_prop[rownames(prelim_anno),],
         annotation_row = prelim_anno,
         annotation_col = azimuth_anno[,"cellType_broad",drop = FALSE],
         annotation_colors =list(cellType_hc = cell_type_colors, cellType_broad = cell_type_colors_broad),
         # main = "Proportion of HC Cluster",
         cluster_rows = FALSE
         # ,
         # cluster_cols = FALSE
         )
dev.off()

png(here(plot_dir, "azimuth_v_prelim_cluster.png"),height = 1200, width = 800)
pheatmap(prelim_anno_prop[rownames(prelim_anno),],
         annotation_row = prelim_anno,
         annotation_col = azimuth_anno[,"cellType_broad",drop = FALSE],
         annotation_colors =list(cellType_hc = cell_type_colors, cellType_broad = cell_type_colors_broad),
         # main = "Proportion of HC Cluster",
         # cluster_rows = FALSE
         # ,
         # cluster_cols = FALSE
)
dev.off()

## Subset maybe split cell types

## Excit_08

heatmap_prelim <- function(cell_type){
  sce_temp <- sce[,sce$cellType_hc == cell_type]
  sce_temp$cellType_hc <- droplevels(sce_temp$cellType_hc)
  sce_temp$prelimCluster <- droplevels(sce_temp$prelimCluster)
  
  n_nuc = table(sce_temp$prelimCluster)
  row_anno = as.matrix(n_nuc)
  colnames(row_anno) <- "n_nuc"
  
  print(row_anno)
  # 22   32   36   65  180  229  247 
  # 356  662 7720 1512 8757 2765 1253 
  
  prelim_anno <- table(sce_temp$prelimCluster, sce_temp$cellType_azimuth)
  prelim_anno_prop <- sweep(prelim_anno,1,rowSums(prelim_anno),`/`)
  
  png(here(plot_dir, paste0("azimuth_v_prelim-",cell_type,".png")),height = 500, width = 800)
  pheatmap(prelim_anno_prop 
           # ,
           # annotation_row = row_anno
           # annotation_colors =list(cellType_hc = cell_type_colors)
           # main = "Proportion of HC Cluster",
           # cluster_rows = FALSE,
           # cluster_cols = FALSE
  )
  dev.off()
}

split_cell_types <- c("Oligo_01","Excit_08","Inhib_01","Excit_02",'Inhib_05')
map(split_cell_types[5], heatmap_prelim)

## Oligo_1 Deep Dive
sce_oligo1 <- sce[,sce$cellType_hc == "Oligo_01"]
sce_oligo1$cellType_hc <- droplevels(sce_oligo1$cellType_hc)
sce_oligo1$prelimCluster <- droplevels(sce_oligo1$prelimCluster)

table(sce_oligo1$prelimCluster)
table(sce_oligo1$prelimCluster, sce_oligo1$Sample)
# 22   32   36   65  180  229  247 
# 356  662 7720 1512 8757 2765 1253 

prelim_anno_Oligo01 <- table(sce_oligo1$prelimCluster, sce_oligo1$cellType_azimuth)
prelim_anno_Oligo01_prop <- sweep(prelim_anno_Oligo01,1,rowSums(prelim_anno_Oligo01),`/`)

n_nuc = table(sce_oligo1$prelimCluster)
row_anno = as.matrix(n_nuc)
colnames(row_anno) <- "n_nuc"

rownames(prelim_anno_Oligo01_prop) <- as.character(rownames(prelim_anno_Oligo01_prop))
rownames(row_anno) <- as.character(rownames(row_anno))


png(here(plot_dir, "azimuth_v_prelim_Oligo01.png"),height = 500, width = 800)
pheatmap(prelim_anno_Oligo01_prop,
         annotation_row = row_anno
         # annotation_colors =list(cellType_hc = cell_type_colors)
         # main = "Proportion of HC Cluster",
         # cluster_rows = FALSE,
         # cluster_cols = FALSE
)
dev.off()

my_plotMarkers(sce = sce_oligo1, 
               marker_list = markers.mathys.tran,
               cat = "prelimCluster",
               pdf_fn = here(plot_dir, "Oligo01_prelim_mathys_markers.pdf"))



## Rand Index
pairwiseRand(sce$cellType_hc, sce$cellType_azimuth, mode="index")
# [1] 0.3094603
pairwiseRand(sce$cellType_k, sce$cellType_azimuth, mode="index")
# [1] 0.2558036
pairwiseRand(sce$cellType_hc, sce$cellType_azimuth, mode="index")
# [1] 0.2441774

library(ggplot2)

## load the query
# p1 <- DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = TRUE) + labs(title = "DLPFC Data")
# # +
# #   NoLegend()
# p2 <- DimPlot(query, group.by = "cellType_hc", label = TRUE, label.size = 3, repel = TRUE) +
#   NoLegend()
# p3 <- DimPlot(query, group.by = "cellType_broad_hc", label = TRUE, label.size = 3, repel = TRUE) +
#   NoLegend()
# p4 <- DimPlot(query, group.by = "Sample")
# 
# ggsave((p1 + p2)/(p3 + p4), filename = here(plot_dir, "Azimuth_UMAP.png"), width = 13, height = 10)
# 
# #### What does the refrence UMAP look like? ####
# reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")
# ref_umap <- DimPlot(reference$map, group.by = "subclass", label = TRUE, label.size = 3, repel = TRUE) +labs(title = "Motor Cortex Reference")
# ggsave(ref_umap + p1, filename = here(plot_dir, "UMAP_Azimuth_reference.png"), width = 14)
# 


#### OUR UMAP space ####
sce$azimuth_cellType <- query$predicted.subclass

library(scater)

UMAP_azi_cellTypes <- ggcells(sce, mapping=aes(x=UMAP.1, y=UMAP.2, colour=azimuth_cellType)) +
  geom_point(size = 0.2, alpha = 0.3) +
  # scale_color_manual(values = cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
  theme_bw() +
  coord_equal()

ggsave(UMAP_azi_cellTypes + 
         guides(colour = guide_legend(override.aes = list(size=2, alpha = 1))) +
         labs(title = "DLPFC UMAP - Azimuth Annotation"),
       filename = here(plot_dir, "UMAP_azi_cellType.png"), width = 10)

ggsave(UMAP_azi_cellTypes + theme(legend.position = "None") + labs(title = "DLPFC UMAP - Azimuth Annotation") + 
         UMAP_azi_cellTypes + theme(legend.position = "None") + facet_wrap(~azimuth_cellType),
       filename = here(plot_dir, "UMAP_azi_cellType_facet.png"), width = 12)


UMAP_hc_cellTypes <- ggcells(sce, mapping=aes(x=UMAP.1, y=UMAP.2, colour=cellType_hc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values =  metadata(sce)$cell_type_colors[levels(sce$cellType_hc)], drop = TRUE) +
  theme_bw() +
  coord_equal()

ggsave(UMAP_hc_cellTypes + theme(legend.position = "None") + labs(title = "DLPFC UMAP - Hierarchical Cluster") + 
         UMAP_hc_cellTypes + theme(legend.position = "None") + facet_wrap(~cellType_hc),
       filename = here(plot_dir, "UMAP_hc_cellType_facet.png"), width = 12)

#### Marker Plots ####

my_plotMarkers(sce = sce, 
               marker_list = markers.mathys.tran,
               cat = "azimuth_basic",
               pdf_fn = here(plot_dir, "Azimuth_basic_mathys_markers.pdf"))


# sgejobs::job_single('05_azimuth_validation', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript 05_azimuth_validation.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
