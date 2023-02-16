library("SingleCellExperiment")
library("scater")
# library("tidyverse")
library("here")
library("sessioninfo")
library("plotly")

## Trying
library("tidySingleCellExperiment")

#### set up plotting ####
source(here("code", "03_build_sce", "my_plotExpression.R"))

plot_dir <- here("plots", "05_explore_sce", "04_3D_UMAP")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

#### Load data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

message("running 3D UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY", ncomponents = 3)
message("done - ", Sys.time())

UMAP_3D <- sce |>
    plot_ly(
        x = ~`UMAP1`,
        y = ~`UMAP2`,
        z = ~`UMAP3`,
        color = ~cellType_hc,
        colors = metadata(sce)$cell_type_colors
    ) |>
    add_markers(size = I(1))

library(htmlwidgets)
saveWidget(UMAP_3D, file = here(plot_dir, "UMAP_3D.html"))
