######### Cell Annotation with Azimuth#####
# devtools::install_github("satijalab/seurat-data") # install SeuratData package
# devtools::install_github("satijalab/azimuth", ref = "release/0.4.5") # install Azimuth package

library("SingleCellExperiment")
library("Seurat")
library("Azimuth")
library("SeuratData")
library("patchwork")
library("here")
library("pheatmap")


#### Plot Setup ####
plot_dir = here("plots","05_explore_sce","05_azimuth_validation")
if(!dir.exists(plot_dir)) dir.create(plot_dir)


## Cell Annotation with RunAzimuth####

## Use Azimuth Refrence ##
## new data
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

query <- CreateSeuratObject(counts = as.matrix(counts(sce)), 
                            meta.data = data.frame(colData(sce)),
                            project = "DLPFC")

query <- RunAzimuth(query, reference = "humancortexref") ## Cell annotation with Azimuth

table(query$predicted.subclass)
# Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT      L6 IT L6 IT Car3        L6b 
# 6744       9888      21440        109       6775        345       1158       1155        332        875 
# Lamp5  Micro-PVM      Oligo        OPC      Pvalb       Sncg        Sst  Sst Chodl        Vip       VLMC 
# 802       2507      11384       1836       1513        243       4689         87       4178       1544 

## Save 
# SeuratDisk::SaveH5Seurat(query, filename = here("processed-data","05_explore_sce","05_azimuth_validation","sce_DLPFC.h5Seurat"))

#### Compare to HC annotations 
ct_annos <- table(query$cellType_hc, query$predicted.subclass)

png(here(plot_dir, "azimuth_v_hc.png"),height = 800, width = 800)
pheatmap(ct_annos,
         main = "Number of Nuclei",
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

## Prop of HC annotation groups
ct_annos_prop <- sweep(ct_annos,1,rowSums(ct_annos),`/`)

ct_max <- apply(ct_annos_prop, 1, max)
table(ct_max  > 0.7)

sum(ct_annos_prop > .90)

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

library(bluster)
jacc.mat <- linkClustersMatrix(query$cellType_hc, query$predicted.subclass)

png(here(plot_dir, "azimuth_v_hc-jacc.png"),height = 800, width = 800)
pheatmap(jacc.mat,
         main = "Strength of the correspondence between Azimuth & HC")
dev.off()

pairwiseRand(query$cellType_hc, query$predicted.subclass, mode="index")
# [1] 0.3094603
pairwiseRand(query$cellType_broad_hc, query$predicted.subclass, mode="index")
# [1] 0.2558036
pairwiseRand(query$cellType_k, query$predicted.subclass, mode="index")
# [1] 0.2441774

library(ggplot2)
p1 <- DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = TRUE) + labs(title = "DLPFC Data")
# +
#   NoLegend()
p2 <- DimPlot(query, group.by = "cellType_hc", label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
p3 <- DimPlot(query, group.by = "cellType_broad_hc", label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
p4 <- DimPlot(query, group.by = "Sample")

ggsave((p1 + p2)/(p3 + p4), filename = here(plot_dir, "Azimuth_UMAP.png"), width = 13, height = 10)

#### What does the refrence UMAP look like? ####
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")
ref_umap <- DimPlot(reference$map, group.by = "subclass", label = TRUE, label.size = 3, repel = TRUE) +labs(title = "Motor Cortex Reference")
ggsave(ref_umap + p1, filename = here(plot_dir, "UMAP_Azimuth_reference.png"), width = 14)



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
load(here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"), verbose = TRUE)


source(here("code", "03_build_sce","my_plotExpression.R"))

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

