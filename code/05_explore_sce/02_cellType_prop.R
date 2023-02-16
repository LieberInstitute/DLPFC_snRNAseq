library("SingleCellExperiment")
library("DeconvoBuddies")
library("tidyverse")
library("patchwork")
library("here")

#### plot setup  ####
my_theme <- theme_bw() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "05_explore_sce", "02_cellType_prop")
# load(here("processed-data", "03_build_sce","cell_type_colors.Rdata"), verbose = TRUE)

#### Load SCE ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## get proportions before dropping ambig
prop_ambig <- as.data.frame(colData(sce)) |>
    group_by(cellType_hc, Sample, Position) |>
    summarize(n = n()) |>
    group_by(Sample) |>
    mutate(prop = n / sum(n))

cell_type_colors_ambig <- metadata(sce)$cell_type_colors[levels(prop_ambig$cellType_hc)]

## Exclude ambig cells
sce <- sce[, sce$cellType_hc != "Ambiguous"]
sce$cellType_hc <- droplevels(sce$cellType_hc)
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)

cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]
cell_type_colors_broad <- metadata(sce)$cell_type_colors[levels(sce$cellType_broad_hc)]
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]

pd <- as.data.frame(colData(sce))

table(pd$cellType_broad_hc)
# Astro EndoMural     Micro     Oligo       OPC     Excit     Inhib
# 3979      2157      1601     10894      1940     24809     11067

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
    labs(y = "Proportion") +
    scale_fill_manual(values = cell_type_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(sample_prop_bar, filename = here(plot_dir, "prop_bar_Sample.png"), width = 12)

## Facet by position
sample_prop_bar_position <- ggplot(data = prop_all, aes(x = Sample, y = prop, fill = cellType_hc)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")),
        size = 2.7,
        position = position_stack(vjust = 0.5)
    ) +
    facet_grid(. ~ Position, scales = "free", space = "free") +
    scale_fill_manual(values = cell_type_colors) +
    labs(y = "Cell Type Proportion", fill = "Cell Type") +
    my_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sample_prop_bar_position, filename = here(plot_dir, "prop_bar_Sample_position.png"), width = 12)
ggsave(sample_prop_bar_position, filename = here(plot_dir, "prop_bar_Sample_position.pdf"), width = 12)

sample_prop_bar_minimal <- ggplot(data = prop_all, aes(x = Sample, y = prop, fill = cellType_hc)) +
    geom_bar(stat = "identity") +
    # facet_grid(. ~ Position, scales = "free", space = "free") +
    scale_fill_manual(values = cell_type_colors) +
    my_theme +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "None"
    )

ggsave(sample_prop_bar_minimal, filename = here(plot_dir, "prop_bar_Sample_minimal.png"), height = 5, width = 6)
ggsave(sample_prop_bar_minimal, filename = here(plot_dir, "prop_bar_Sample_minimal.pdf"), height = 5, width = 6)


## number nuclei
sample_n_bar_position <- ggplot(data = prop_all, aes(x = Sample, y = n, fill = cellType_hc)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(prop > 0.02, format(round(n, 3), 3), "")),
        size = 2.7,
        position = position_stack(vjust = 0.5)
    ) +
    facet_grid(. ~ Position, scales = "free", space = "free") +
    scale_fill_manual(values = cell_type_colors) +
    labs(y = "n Nuclei", fill = "Cell Type") +
    my_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sample_n_bar_position, filename = here(plot_dir, "nnuc_bar_Sample_position.png"), width = 12)
ggsave(sample_n_bar_position, filename = here(plot_dir, "nnuc_bar_Sample_position.pdf"), width = 12)


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

#### Plots ambig Plots ####
prop_ambig_plus <- prop_ambig |>
    mutate(ambig = "Pre-drop") |>
    bind_rows(prop_all |> mutate(ambig = "Post-drop"))


prop_ambig_bar <- ggplot(data = prop_ambig_plus, aes(x = Sample, y = prop, fill = cellType_hc)) +
    geom_bar(stat = "identity") +
    geom_text(
        aes(
            label = ifelse(prop > 0.02, format(round(prop, 3), 3), ""),
            color = as.character(cellType_hc) == "ambig"
        ),
        # geom_text(aes(label = ifelse(prop > 0.02, cellType, "")),
        size = 2.5,
        position = position_stack(vjust = 0.5)
        # color = "gray35"
    ) +
    scale_fill_manual(values = c(cell_type_colors_ambig)) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    labs(y = "Proportion", fill = "Cell Type") +
    facet_grid(fct_rev(ambig) ~ Position, scales = "free", space = "free") +
    # my_theme +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(color = FALSE, fill = guide_legend(ncol = 1))

ggsave(prop_ambig_bar, filename = here(plot_dir, "prop_bar_ambig_Position.png"), width = 11, height = 9)
ggsave(prop_ambig_bar, filename = here(plot_dir, "prop_bar_ambig_Position.pdf"), width = 11, height = 9)

#### Layer Annotation Proportions ####
table(pd$cellType_layer)
# Astro    EndoMural        Micro        Oligo          OPC   Excit_L2/3     Excit_L3 Excit_L3/4/5     Excit_L4
# 3979         2157         1601        32051         1940           82        10459         3043         2388
# Excit_L5   Excit_L5/6     Excit_L6        Inhib
# 2505         2487         1792        11067
table(pd$layer_annotation)

n_nuc_layer <- pd |>
    group_by(cellType_hc, cellType_layer) |>
    summarize(n_nuc = n()) |>
    mutate(
        ct_short = sub("^(\\w).*(_\\d)", "\\1\\2", cellType_hc),
        anno = case_when(
            ct_short == cellType_hc ~ as.character(n_nuc),
            n_nuc < 2000 ~ paste(ct_short, n_nuc),
            TRUE ~ paste0(ct_short, "\n", n_nuc)
        )
    )

## Overlap between
# Excit_14 82
# Excit_15 66
## Solve with ifelse..add 66 in post?
barplot_n_nuc_layer <- n_nuc_layer |>
    ggplot(aes(x = cellType_layer, y = n_nuc, fill = cellType_hc)) +
    geom_bar(stat = "identity") +
    # geom_text(aes(label = ifelse(n_nuc > 70, n_nuc, NA)), size = 2.5, position = position_stack(vjust = 0.5)) +
    geom_text(aes(label = anno), size = 2.5, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = cell_type_colors) +
    my_theme +
    theme(legend.position = "None", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Number of Nuclei", x = "Layer Cell Type Annotation")

ggsave(barplot_n_nuc_layer, filename = here(plot_dir, "barplot_n_nuc_layer.png"), width = 8)
ggsave(barplot_n_nuc_layer, filename = here(plot_dir, "barplot_n_nuc_layer.pdf"), width = 8)

## Prop Bar
prop_layer <- pd |>
    # filter(!is.na(cellType_layer)) |>
    count(Sample, cellType_layer) |>
    left_join(n_nuc) |>
    mutate(
        prop = n / n_nuc,
        cellType_layer = addNA(cellType_layer)
    )

layer_prop_bar_all <- plot_composition_bar(prop_layer,
    sample_col = "Sample",
    ct_col = "cellType_layer",
    min_prop_text = .02
) +
    scale_fill_manual(values = cell_type_colors_layer) +
    theme_bw()
ggsave(layer_prop_bar_all, filename = here(plot_dir, "prop_bar_layer_ALL.png"), width = 3)

#### compare hc v layer ####
prop_compare <- pd |>
    group_by(cellType_layer) |>
    summarize(prop = n() / nrow(pd)) |>
    mutate(
        Annotation = "Layer",
        cellType = as.character(cellType_layer)
    ) |>
    bind_rows(pd |>
        group_by(cellType_hc, cellType_layer) |>
        summarize(prop = n() / nrow(pd)) |>
        mutate(
            Annotation = "hc",
            cellType = as.character(cellType_hc)
        )) |>
    replace_na(list(cellType = "Exclude")) |>
    group_by(cellType_layer) |>
    mutate(
        order = as.numeric(cellType_layer),
        order2 = row_number()
    ) |>
    replace_na(list(order = 12.5)) |> ## cheat Excluded cells to be between Excit & Inhib
    ungroup() |>
    mutate(
        order = order + (order2 * .1),
        cellType = factor(cellType, levels = unique(cellType[order(order)]))
    )

prop_compare |>
    filter(Annotation == "hc") |>
    arrange(cellType_layer) |>
    print(n = 29)
prop_compare |> filter(is.na(cellType_layer))


prop_compare_bar <- ggplot(data = prop_compare, aes(x = Annotation, y = prop, fill = cellType)) +
    geom_bar(stat = "identity", width = .99) +
    geom_text(aes(label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")),
        # geom_text(aes(label = ifelse(prop > 0.02, cellType, "")),
        size = 3,
        position = position_stack(vjust = 0.5)
    ) +
    scale_fill_manual(values = c(cell_type_colors, cell_type_colors_layer)) +
    labs(y = "Proportion", fill = "Cell Type") +
    my_theme +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.01, 0.01))) + ## reduce white space in boarder
    theme(
        panel.grid = element_blank(),
        # panel.border = element_blank(),
        legend.position = "None"
    )

ggsave(prop_compare_bar, filename = here(plot_dir, "prop_compare_bar.png"), width = 2, height = 8)
ggsave(prop_compare_bar, filename = here(plot_dir, "prop_compare_bar.pdf"), width = 2, height = 8)


# prop_bar_fineVlayer <- prop_bar_all + theme(legend.position = "None") + labs(y = "Mean Proption - hc annotation k=29") +
#     layer_prop_bar_all + labs(y = "Mean Proption - layer annotation k=13")
#
# ggsave(prop_bar_fineVlayer, filename = here(plot_dir, "prop_bar_fineVlayer_ALL.png"))
# ggsave(prop_bar_fineVlayer, filename = here(plot_dir, "prop_bar_fineVlayer_ALL.pdf"))

## Sample prop bar
layer_sample_prop_bar <- plot_composition_bar(prop_layer,
    sample_col = "Sample",
    x_col = "Sample",
    ct_col = "cellType_layer",
    min_prop_text = .03
) +
    labs(y = "Proportion") +
    scale_fill_manual(values = metadata(sce)$cell_type_colors_layer) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(layer_sample_prop_bar, filename = here(plot_dir, "prop_bar_layer_Sample.png"), width = 12)
ggsave(layer_sample_prop_bar, filename = here(plot_dir, "prop_bar_layer_Sample.pdf"), width = 12)

## facet version
layer_prop_bar_position <- ggplot(data = prop_layer, aes(x = Sample, y = prop, fill = cellType_layer)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = ifelse(prop > 0.02, format(round(prop, 3), 3), "")),
        size = 2.7,
        position = position_stack(vjust = 0.5)
    ) +
    facet_grid(. ~ Position, scales = "free", space = "free") +
    scale_fill_manual(values = metadata(sce)$cell_type_colors_layer) +
    labs(y = "Cell Type Proportion", fill = "Cell Type") +
    my_theme +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None"
    )

ggsave(layer_prop_bar_position, filename = here(plot_dir, "prop_bar_layer_position.png"), width = 10)
ggsave(layer_prop_bar_position, filename = here(plot_dir, "prop_bar_layer_position.pdf"), width = 10)



# sgejobs::job_single('02_cellType_props', create_shell = TRUE, queue= 'bluejay', memory = '10G', command = "Rscript 02_cellType_prop.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
