library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("patchwork")
library("here")

#### Load Data ####
# load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
sce <- HDF5Array::loadHDF5SummarizedExperiment(here("processed-data", "sce", "sce_DLPFC_annotated"))

#### Set Up Plotting ####
# my_theme <- theme_classic() +
#     theme(text = element_text(size = 15),
#           strip.background = element_rect(colour = "black", fill = "light grey"))

my_theme <- theme_bw() +
    theme(
        text = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )


plot_dir <- here("plots", "05_explore_sce", "01_reduced_dim_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### TSNE with drop nuc ####
levels(sce$cellType_hc)
sce$cellType_hc <- forcats::fct_relevel(sce$cellType_hc, "drop", after = Inf)

cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]

TSNE_cellTypes_hc <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cell_type_colors) +
  my_theme +
  coord_equal() +
  labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")


## No legends
TSNE_cellTypes_hc_no_legend <- TSNE_cellTypes_hc + theme(legend.position = "None")

TSNE_cellTypes_hc_facet <- TSNE_cellTypes_hc_no_legend +
  facet_wrap(~cellType_hc)

## Add full + facet for cell types
ggsave(
  TSNE_cellTypes_hc_no_legend +
    TSNE_cellTypes_hc_facet +
    theme(axis.title.y = element_blank()),
  filename = here(plot_dir, "TSNE_cellType_full_facet-drop.png"),
  width = 13
)

ggsave(
  TSNE_cellTypes_hc_no_legend +
    TSNE_cellTypes_hc_facet +
    theme(axis.title.y = element_blank()),
  filename = here(plot_dir, "TSNE_cellType_full_facet-drop.pdf"),
  width = 13
)

## Exclude drop nuc
sce <- sce[, sce$cellType_hc != "drop"]
sce$cellType_hc <- droplevels(sce$cellType_hc)

## Adjust color pallets
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
ggsave(
    UMAP_cellTypes_hc +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "UMAP_cellType.png"), width = 9
)
ggsave(
    UMAP_cellTypes_hc +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "UMAP_cellType.pdf"), width = 9
)

## No legends
UMAP_cellTypes_hc_no_legend <- UMAP_cellTypes_hc + theme(legend.position = "None")

ggsave(UMAP_cellTypes_hc_no_legend,
    filename = here(plot_dir, "UMAP_cellType_no_legend.png")
)
ggsave(UMAP_cellTypes_hc_no_legend,
    filename = here(plot_dir, "UMAP_cellType_no_legend.pdf")
)

## Add facet for cell types
UMAP_cellTypes_hc_facet <- UMAP_cellTypes_hc_no_legend +
    facet_wrap(~cellType_hc) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3)

ggsave(UMAP_cellTypes_hc_facet,
    filename = here(plot_dir, "UMAP_cellType_facet.png"),
)
ggsave(UMAP_cellTypes_hc_facet,
    filename = here(plot_dir, "UMAP_cellType_facet.pdf"),
)

## Add full + facet for cell types
ggsave(
    UMAP_cellTypes_hc_no_legend +
        UMAP_cellTypes_hc_facet +
        theme(axis.title.y = element_blank()),
    filename = here(plot_dir, "UMAP_cellType_full_facet.png"),
    width = 13
)

ggsave(
    UMAP_cellTypes_hc_no_legend +
        UMAP_cellTypes_hc_facet +
        theme(axis.title.y = element_blank()),
    filename = here(plot_dir, "UMAP_cellType_full_facet.pdf"),
    width = 13
)

#### Plot clusters in TSNE ####
TSNE_cellTypes_hc <- ggcells(sce, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_hc)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors) +
    my_theme +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

## Add large legends
ggsave(
    TSNE_cellTypes_hc +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_cellType.png"), width = 9
)
ggsave(
    TSNE_cellTypes_hc +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_cellType.pdf"), width = 9
)

## No legends
TSNE_cellTypes_hc_no_legend <- TSNE_cellTypes_hc + theme(legend.position = "None")

ggsave(TSNE_cellTypes_hc_no_legend,
    filename = here(plot_dir, "TSNE_cellType_no_legend.png")
)
ggsave(TSNE_cellTypes_hc_no_legend,
    filename = here(plot_dir, "TSNE_cellType_no_legend.pdf")
)

## Add facet for cell types
TSNE_cellTypes_hc_facet <- TSNE_cellTypes_hc_no_legend +
    facet_wrap(~cellType_hc)

ggsave(TSNE_cellTypes_hc_facet,
    filename = here(plot_dir, "TSNE_cellType_facet.png"),
)
ggsave(TSNE_cellTypes_hc_facet,
    filename = here(plot_dir, "TSNE_cellType_facet.pdf"),
)

## Add full + facet for cell types
ggsave(
    TSNE_cellTypes_hc_no_legend +
        TSNE_cellTypes_hc_facet +
        theme(axis.title.y = element_blank()),
    filename = here(plot_dir, "TSNE_cellType_full_facet.png"),
    width = 13
)

ggsave(
    TSNE_cellTypes_hc_no_legend +
        TSNE_cellTypes_hc_facet +
        theme(axis.title.y = element_blank()),
    filename = here(plot_dir, "TSNE_cellType_full_facet.pdf"),
    width = 13
)

#### TSNE layer annotations ####
TSNE_cellType_layer <- ggcells(sce[, !is.na(sce$cellType_layer)], mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType_layer)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_type_colors_layer) +
    my_theme +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2")

ggsave(
    TSNE_cellType_layer +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_cellType_layer.png"), width = 9
)

ggsave(
    TSNE_cellType_layer +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))),
    filename = here(plot_dir, "TSNE_cellType_layer.pdf"), width = 9
)

## no legend
ggsave(
    TSNE_cellType_layer +
        theme(legend.position = "None"),
    filename = here(plot_dir, "TSNE_cellType_layer_no_legend.png")
)
ggsave(
    TSNE_cellType_layer +
        theme(legend.position = "None"),
    filename = here(plot_dir, "TSNE_cellType_layer_no_legend.pdf")
)



# sgejobs::job_single('01_reduce_dim_plots', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 01_reduce_dim_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
