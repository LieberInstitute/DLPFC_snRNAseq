library(liana)
library(tidyverse)
library(unglue)
library(SingleCellExperiment)
library(here)
library(ComplexHeatmap)
library(circlize)
library(glue)

# sce color ---------------------------------------------------------------
sce <- readRDS(here("processed-data/sce/sce_DLPFC_annotated/se.rds"))
cell_type_colors_layer <- metadata(sce)$cell_type_colors_layer[levels(sce$cellType_layer)]
# cell_type_colors_broad <- metadata(sce)$cell_type_colors_broad[levels(sce$cellType_broad_hc)]


# factor_cell_type_broad <- function(vec){
#     factor(vec,
#            levels = c("astro",  "endomural", "micro", "oligo", "opc", "excit", "inhib"
#            ),
#            labels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC", "Excit", "Inhib" )
#     )
# }
#
# factor_cell_type_layer <- function(vec){
#     factor(vec,
#            # levels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
#            #            "Excit_L2_3", "Excit_L3", "Excit_L3_4_5", "Excit_L4","Excit_L5",
#            #            "Excit_L5_6","Excit_L6","Inhib"),
#            levels = c("astro",  "endomural", "micro", "oligo", "opc",
#                       "excit_l2_3", "excit_l3", "excit_l3_4_5", "excit_l4","excit_l5",
#                       "excit_l5_6","excit_l6","inhib"),
#            labels = c("Astro",  "EndoMural", "Micro", "Oligo", "OPC",
#                       "Excit_L2/3", "Excit_L3", "Excit_L3/4/5", "Excit_L4","Excit_L5",
#                       "Excit_L5/6","Excit_L6","Inhib")
#
#     )
# }


# NOTE: this line only relevents to Boyi
# .libPaths(c(.libPaths(), "/users/bguo/R/4.2"))

plot_path <- here("plots/07_ccc/")

# Per-Section lots -------------------------------------------------------
scn_fld_path <- here("processed-data/07_ccc/")
# crn_scn <- c("Anterior", "Middle", "Posterior")
# crn_scn |>
# walk(.f = function(scn){
liana_res <- readRDS(here(
    scn_fld_path,
    "liana_consensus.rds"
))
liana_trunc <- liana_res %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01)


liana_trunc |>
    arrange(aggregate_rank) |>
    write.csv(
        file = here(
            "processed-data/07_ccc/",
            "TableS9_LIANA_LR.csv"
        ),
        row.names = FALSE
    )

nrow(liana_trunc)
# 6588

liana_trunc |>
    transmute(
        pair = glue::glue("{ligand.complex}_{receptor.complex}")
    ) |>
    summarize(n_distinct(pair))


# Create Folder
if (!dir.exists(plot_path)) {
    dir.create(plot_path, recursive = TRUE)
}

# Save Heatmap
pdf(
    file = here(plot_path, "/liana_heatmap_overall.pdf"),
    width = 9,
    height = 8
)
plot_liana_heatmap(
    liana_trunc |>
        mutate(
            source = factor(
                source,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            ),
            target = factor(
                target,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            )
        ) |>
        liana:::.get_freq(),
    cell_col = cell_type_colors_layer # ,
    # row_dend_reorder = FALSE,
    # column_dend_reorder = FALSE
) |>
    print()
dev.off()
# })


# Circular plot ------------------------------------------------------------


lig_gene <- "EFNA5"
rec_gene <- "EPHA5"

lig_rec_df <- liana_trunc |>
    filter(
        ligand.complex == lig_gene,
        receptor.complex == rec_gene
    )


full_join(
    table(lig_rec_df$source) |> prop.table() |> data.frame(),
    table(lig_rec_df$target) |> prop.table() |> data.frame(),
    by = "Var1",
    suffix = c("_source", "_target")
)

# xtabs(~source + target, data = tmp)


# pdf(filename = here(plot_path, "/liana_circ_plot_EFNA5-EPHA5.pdf")#,
#     # width = 9,
#     # height = 8
#     )
png(
    file = here(plot_path, "/liana_circ_plot_EFNA5-EPHA5.png"),
    height = 1200,
    width = 1200, units = "px"
)

lig_rec_df |>
    mutate(
        source = factor(
            source,
            levels = c(
                "Astro", "EndoMural", "Micro",
                "Oligo", "OPC", "Excit_L2/3",
                "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
            )
        ),
        target = factor(
            target,
            levels = c(
                "Astro", "EndoMural", "Micro",
                "Oligo", "OPC", "Excit_L2/3",
                "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
            )
        )
    ) |>
    plot_liana_circ(
        cell_col = cell_type_colors_layer
    ) |>
    print()

dev.off()

# tmp |> liana_dotplot()

# mutate(LR_pair = glue::glue("{ligand.complex}-{receptor.complex}")) |>
#     pull(LR_pair) |>
#     unique()



# Plot Pie Chart for Sender -----------------------------------------------

ggsave(
    filename = here(plot_path, "liana_pie_EFNA5-EPHA5_sender_nolegend.pdf"),
    plot = lig_rec_df |>
        mutate(
            source = factor(
                source,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            ),
            target = factor(
                target,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            )
        ) |>
        ggplot() +
        geom_bar(aes(x = factor(1), fill = source),
            stat = "count",
            show.legend = FALSE
        ) +
        coord_polar("y", start = 0) +
        scale_fill_manual(
            values = cell_type_colors_layer
        ) +
        labs(title = "Source") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 20))
    # theme_bw() +
    #     theme(axis.text.y = element_blank(),
    #           axis.title.y = element_blank(),
    #           axis.ticks.y = element_blank(),
    #           axis.title.x = element_blank(),
    #           axis.text.x = element_blank())
)

# Plot Pie Chart for Target -----------------------------------------------

ggsave(
    filename = here(plot_path, "liana_pie_EFNA5-EPHA5_target.pdf"),
    plot = lig_rec_df |>
        mutate(
            source = factor(
                source,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            ),
            target = factor(
                target,
                levels = c(
                    "Astro", "EndoMural", "Micro",
                    "Oligo", "OPC", "Excit_L2/3",
                    "Excit_L3", "Excit_L3/4/5", "Excit_L4",
                    "Excit_L5", "Excit_L5/6", "Excit_L6", "Inhib"
                )
            )
        ) |>
        ggplot() +
        geom_bar(aes(x = factor(1), fill = target),
            stat = "count"
            # , show.legend = FALSE
        ) +
        coord_polar("y", start = 0) +
        scale_fill_manual(
            values = cell_type_colors_layer
        ) +
        labs(title = "Target") +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 20))
    # theme_bw() +
    #     theme(axis.text.y = element_blank(),
    #           axis.title.y = element_blank(),
    #           axis.ticks.y = element_blank(),
    #           axis.title.x = element_blank(),
    #           axis.text.x = element_blank())
)


heat_freq()

# Per-Sample lots -------------------------------------------------------
# smp_fld_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/07_ccc/CCC_snRNA_R3.1/"
# file_dic <- data.frame(fl_name = list.files(smp_fld_path)) |>
#     unglue_unnest(fl_name, "{Br}_{sec}", remove = FALSE) |>
#     filter(substr(fl_name, start = 1, stop = 2) == "Br") |>
#     left_join(
#         y = data.frame(
#             sec = c("ant", "mid", "post"),
#             Section = c("Anterior", "Middle", "Posterior")
#         )
#     ) |>
#     select(fl_name, Section)
#
# walk2(.x = file_dic$fl_name,
#       .y = file_dic$Section,
#       .f = function(name, scn){
#           liana_res <- readRDS(paste0(smp_fld_path, name, "/liana_consensus.rds"))
#           liana_trunc <- liana_res %>%
#               # only keep interactions concordant between methods
#               filter(aggregate_rank <= 0.01)
#
#           # Create Folder
#           if(!dir.exists(paste0(plot_path, scn)))
#               dir.create(paste0(plot_path,scn), recursive = TRUE)
#
#           # Save Heatmap
#           jpeg(filename =  paste0(plot_path, scn, "/",name,".jpeg"))
#           print(heat_freq(liana_trunc))
#           dev.off()
#       }
# )
