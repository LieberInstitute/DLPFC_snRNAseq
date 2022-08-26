
library("SingleCellExperiment")
library("DeconvoBuddies")
library("tidyverse")
library("here")

#### plot setup  ####
plot_dir <- here("plots", "05_explore_sce", "02_cellType_prop")
# load(here("processed-data", "03_build_sce","cell_type_colors.Rdata"), verbose = TRUE)
#


#### Load SCE ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]

pd <- as.data.frame(colData(sce))

table(pd$cellType_broad_hc)
# Astro Endo.Mural      Excit      Inhib      Micro      Oligo        OPC
# 3979       2157      24809      11067       1601      32051       1940

n_nuc <- pd |>
    group_by(Sample, BrNum, round, Position, sex, age) |>
    summarize(n_nuc = n()) |>
    ungroup()

prop_all <- pd |>
    count(Sample, cellType_hc) |>
    left_join(n_nuc) |>
    mutate(prop = n / n_nuc)

prop_broad <- pd |>
    count(Sample, cellType_broad_hc) |>
    left_join(n_nuc) |>
    mutate(prop = n / n_nuc)

## n nuc bar plot
barplot_n_nuc <- pd |>
    group_by(cellType_hc) |>
    summarize(n_nuc = n()) |>
    ggplot(aes(x = cellType_hc, y = n_nuc, fill = cellType_hc)) +
    geom_col() +
    geom_text(aes(label = n_nuc), size = 2.5) +
    scale_fill_manual(values = cell_type_colors) +
    theme_bw() +
    theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank()) +
    labs(y = "Number of Nuclei")

ggsave(barplot_n_nuc, filename = here(plot_dir, "barplot_n_nuc.png"), heigh = 3.5, width = 8)



prop_boxplots <- prop_broad |>
    ggplot(aes(x = cellType_broad_hc, y = prop, fill = cellType_broad_hc)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, colour = "black", pch = 21) +
    scale_fill_manual(values = cell_type_colors_broad) +
    # scale_color_manual(values = cell_type_colors_broad) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(prop_boxplots, filename = here(plot_dir, "prop_boxplots.png"))

prop_boxplot_position <- prop_broad |>
    separate(Sample, into = c("BrNum", "Position")) |>
    ggplot(aes(x = Position, y = prop, fill = cellType_broad_hc)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.2, colour = "black", pch = 21) +
    scale_fill_manual(values = cell_type_colors_broad) +
    facet_wrap(~cellType_broad_hc, scales = "free_y") +
    theme_bw() +
    theme(legend.position = "None")

ggsave(prop_boxplot_position, filename = here(plot_dir, "prop_boxplot_position.png"))


#### Compositon Plots ####
library(patchwork)
## Mean Prop
prop_bar_all <- plot_composition_bar(prop_all,
    sample_col = "Sample",
    ct_col = "cellType_hc",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors) +
    theme_bw()

ggsave(prop_bar_all, filename = here(plot_dir, "prop_bar_ALL.png"))

broad_prop_bar_all <- plot_composition_bar(prop_broad,
    sample_col = "Sample",
    ct_col = "cellType_broad_hc",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors_broad) +
    theme_bw() +
    theme(legend.position = "None")

ggsave(broad_prop_bar_all, filename = here(plot_dir, "prop_bar_broad_ALL.png"))

ggsave(broad_prop_bar_all + theme(legend.position = "None") + labs(y = "Mean Proption - Broad") + prop_bar_all,
    filename = here(plot_dir, "prop_bar_fineVbroad_ALL.png")
)

## By Sample
sample_prop_bar <- plot_composition_bar(prop_all,
    sample_col = "Sample",
    x_col = "Sample",
    ct_col = "cellType_hc",
    min_prop_text = .03
) +
    scale_fill_manual(values = cell_type_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(sample_prop_bar, filename = here(plot_dir, "prop_bar_Sample.png"), width = 12)

broad_prop_bar <- plot_composition_bar(prop_broad,
    sample_col = "Sample",
    x_col = "Sample",
    ct_col = "cellType_broad_hc",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors_broad) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(broad_prop_bar, filename = here(plot_dir, "prop_bar_broad_Sample.png"), width = 12)


## Examine Round
broad_prop_bar_round <- plot_composition_bar(prop_broad,
    sample_col = "Sample",
    x_col = "round",
    ct_col = "cellType_broad_hc",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors_broad)

ggsave(broad_prop_bar_round, filename = here(plot_dir, "prop_bar_broad_round.png"))

## Position/Position
broad_prop_bar_pos <- plot_composition_bar(prop_broad,
    sample_col = "Sample",
    x_col = "Position",
    ct_col = "cellType_broad_hc",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors_broad) +
    theme_bw() +
    labs(x = "Position")

ggsave(broad_prop_bar_pos, filename = here(plot_dir, "prop_bar_broad_Position.png"))
