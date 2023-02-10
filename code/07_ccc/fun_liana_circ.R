plot_liana_circ <- function(liana_res, source_groups = NULL, target_groups = NULL,
    cex = 1, transparency = 0.4,
    facing = "clockwise", adj = c(-0.5, 0.05),
    cell_col = NULL, ...) {
    freqs <- liana_res %>%
        # (if (!is.null(source_groups))
        # filter(., source %in% source_groups)
        # else .) %>% (if (!is.null(target_groups))
        #     filter(., target %in% target_groups)
        #     else .) %>%
        liana:::.get_freq()
    celltypes <- union(colnames(freqs), rownames(freqs))
    if (is.null(cell_col)) {
        grid.col <- (grDevices::colorRampPalette(
            (RColorBrewer::brewer.pal(n = 8, name = "Dark2"))
        )
        )(length(celltypes)) %>%
            setNames(celltypes)
    } else {
        grid.col <- cell_col
    }


    circlize::circos.clear()
    circlize::chordDiagram(freqs,
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        link.arr.type = "big.arrow",
        transparency = transparency,
        grid.col = grid.col,
        annotationTrack = c("grid"), self.link = 1,
        big.gap = 7.5, small.gap = 5, ...
    )
    circlize::circos.trackPlotRegion(
        track.index = 1,
        panel.fun = function(x, y) {
            xlim <- circlize::get.cell.meta.data("xlim")
            ylim <- circlize::get.cell.meta.data("ylim")
            sector.name <- circlize::get.cell.meta.data("sector.index")
            circlize::circos.text(mean(xlim), ylim[1], sector.name,
                facing = facing, niceFacing = TRUE, adj = adj, cex = cex
            )
        }, bg.border = NA
    )
    p <- grDevices::recordPlot()
    return(p)
}
