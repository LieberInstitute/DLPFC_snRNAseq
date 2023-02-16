#' Plot list of marker genes
#'
#' @param sce
#' @param marker_list
#' @param assay
#' @param cat
#' @param fill_colors
#' @param pdf_fn
#'
#' @return
#' @export
#'
#' @examples
my_plotMarkers <- function(sce, marker_list, assay = "logcounts", cat = "cellType", fill_colors = NULL, pdf_fn) {
    message("plotting: ", pdf_fn)
    pdf(pdf_fn, height = 6, width = 8)

    for (i in 1:length(marker_list)) {
        message(names(marker_list)[[i]])

        m <- marker_list[[i]]
        markers <- m[m %in% rownames(sce)]

        # if(m != markers) message("Missing...",paste(m[!m %in% markers], collapse = ", "))
        if (!identical(m, markers)) message("Missing markers...")
        print(
            custom_plotExpression(sce,
                genes = markers,
                title = names(marker_list)[[i]],
                cat = cat,
                fill_colors = fill_colors
            )
        )
    }

    dev.off()
}
