plot_liana_heatmap <- function (mat, font_size = 12, grid_text = FALSE,
                                name = "Frequency",
                                pallette = c("white", "violetred2"),
                                row_title = "Sender (Cell types)",
                                column_title = "Receiver (Cell types)",
                                cell_col = NULL, # Boyi Added parameter
                                ...)
{
    if (grid_text) {
        grid_text <- function(j, i, x, y, width, height, fill) {
            grid_text <- grid.text(sprintf("%d", mat[i, j]),
                                   x, y,
                                   gp = gpar(fontsize = font_size * 0.83))
        }
    }
    else {
        grid_text <- NULL
    }

    # browser()
    cell_anno <- unique(rownames(mat))
    if(is.null(cell_col)){
        fun_color <- grDevices::colorRampPalette(
            RColorBrewer::brewer.pal(n = 8, name = "Dark2")

        )
        cell_anno <- fun_color(length(cell_anno)) %>%
            setNames(cell_anno)
    } else {
        cell_anno <- cell_col
    }

    ha_opts <- list(show_legend = FALSE,
                    show_annotation_name = FALSE,
                    col = list(anno = cell_anno),
                    simple_anno_size = grid::unit(0.25, "cm"))

    column_ha <- exec("HeatmapAnnotation", anno = names(cell_anno),
                      !!!ha_opts)

    row_ha <- exec("rowAnnotation", anno = names(cell_anno),
                   !!!ha_opts)

    column_bar <- ComplexHeatmap::HeatmapAnnotation(
        bar = liana:::.anno_barplot(colSums(mat),
                                    cell_anno,
                                    axis.font.size = font_size * 0.4),
        annotation_name_gp = gpar(fontsize = font_size * 0.5),
        show_legend = FALSE, show_annotation_name = FALSE)

    row_bar <- ComplexHeatmap::rowAnnotation(
        bar2 = liana:::.anno_barplot(rowSums(mat),
                                     cell_anno, font_size * 0.4),
        gp = gpar(fill = cell_anno, col = cell_anno),
        show_legend = FALSE, show_annotation_name = FALSE#,
        # gap = unit(10,"points")
        )

    ComplexHeatmap::Heatmap(
        mat, col = colorRampPalette(pallette)(10),
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_names_side = "left",
        top_annotation = column_bar,
        bottom_annotation = column_ha,
        right_annotation = row_bar,
        left_annotation = row_ha,
        row_title = row_title,
        row_names_gp = gpar(fontsize = font_size),
        row_title_gp = gpar(fontsize = font_size * 1.2),
        column_names_gp = gpar(fontsize = font_size),
        column_title = column_title,
        column_title_gp = gpar(fontsize = font_size * 1.2),
        column_title_side = "bottom",
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = font_size *
                                0.9, fontface = "bold"),
            border = NA, labels_gp = gpar(fontsize = font_size * 0.9),
            grid_width = unit(2, "mm")
        ),
        name = name, cell_fun = grid_text,
        ...)
}
