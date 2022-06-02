
library("patchwork")
#' Plot Reduced Dim of SCE
#' wrapper function for ggcells adding facet wrap to selected category
#' 
#' @param type Type of reduced dimension i.e. TSNE or UMAP
#' @param facet_by Feature to facet and color by, column in colData
#'
#' @return Paired reduced dim and faceted reduced dim plot
#' @export
#'
#' @examples
plot_reducedDim_facet <- function(sce, type = "TSNE", facet_by = "round", title =""){
  
  # main_plot <- ggcells(sce, mapping=aes(x=paste0(type,".1"), y=TSNE.2, colour=round)) +
  
  my_x = paste0(type, ".1")
  my_y = paste0(type, ".2")
  
  main_plot <- ggcells(sce, mapping=aes_string(x=my_x, y=my_y, colour=facet_by)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    my_theme + 
    theme(legend.position = "none")
  
  facet_plot <- main_plot + facet_wrap(as.formula(paste("~", facet_by)))
    
  return(main_plot + labs(title = title) + facet_plot)  
  
}


plot_reducedDim_qc <- function(sce, type = "TSNE", color_by = "sum", title =""){
  
  # main_plot <- ggcells(sce, mapping=aes(x=paste0(type,".1"), y=TSNE.2, colour=round)) +
  
  my_x = paste0(type, ".1")
  my_y = paste0(type, ".2")
  
  main_plot <- ggcells(sce, mapping=aes_string(x=my_x, y=my_y, colour=color_by)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    viridis::scale_color_viridis() + 
    my_theme 
  
  return(main_plot + labs(title = title))  
  
}
