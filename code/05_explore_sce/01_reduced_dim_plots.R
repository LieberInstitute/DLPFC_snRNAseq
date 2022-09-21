library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("here")

#### Load Data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

#### Set Up Plotting ####
my_theme <- theme_classic() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "05_explore_sce", "01_reduced_dim_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]

#### UMAP plots ####
## Cell Type UMAP
UMAP_cellTypes_hc <- ggcells(sce, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors) +
    my_theme +
    coord_equal() +
  labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2")

## Add large legends
ggsave(UMAP_cellTypes_hc +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_cellType.png"), width = 9
)
ggsave(UMAP_cellTypes_hc +
    guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
filename = here(plot_dir, "UMAP_cellType.pdf"), width = 9
)

## No legends
ggsave(UMAP_cellTypes_hc +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_no_legend.png")
)
ggsave(UMAP_cellTypes_hc +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_no_legend.pdf")
)

## Add facet for cell types
ggsave(UMAP_cellTypes_hc +
    facet_wrap(~cellType_hc) +
      scale_x_continuous(n.breaks = 3) +
      scale_y_continuous(n.breaks = 3) +
    theme(legend.position = "None"),
filename = here(plot_dir, "UMAP_cellType_facet.png"), 
)
ggsave(UMAP_cellTypes_hc +
         facet_wrap(~cellType_hc) +
         scale_x_continuous(n.breaks = 3) +
         scale_y_continuous(n.breaks = 3) +
         theme(legend.position = "None"),
       filename = here(plot_dir, "UMAP_cellType_facet.pdf"), 
)


#### Plot clusters in TSNE ####
## Cell Type UMAP
TSNE_cellTypes_hc <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors) +
  my_theme +
  coord_equal() +
  labs(x = "tSNE Dimension 1", y = "tSNE Dimension 2")

## Add large legends
ggsave(TSNE_cellTypes_hc +
         guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
       filename = here(plot_dir, "TSNE_cellType.png"), width = 9
)
ggsave(TSNE_cellTypes_hc +
         guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
       filename = here(plot_dir, "TSNE_cellType.pdf"), width = 9
)

## No legends
ggsave(TSNE_cellTypes_hc +
         theme(legend.position = "None"),
       filename = here(plot_dir, "TSNE_cellType_no_legend.png")
)
ggsave(TSNE_cellTypes_hc +
         theme(legend.position = "None"),
       filename = here(plot_dir, "TSNE_cellType_no_legend.pdf")
)

## Add facet for cell types
ggsave(TSNE_cellTypes_hc +
         facet_wrap(~cellType_hc) +
         # scale_x_continuous(n.breaks = 3) +
         # scale_y_continuous(n.breaks = 3) +
         theme(legend.position = "None"),
       filename = here(plot_dir, "TSNE_cellType_facet.png"), 
)
ggsave(TSNE_cellTypes_hc +
         facet_wrap(~cellType_hc) +
         # scale_x_continuous(n.breaks = 3) +
         # scale_y_continuous(n.breaks = 3) +
         theme(legend.position = "None"),
       filename = here(plot_dir, "TSNE_cellType_facet.pdf"), 
)


#### layer annotaions ####
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

