
plot_reducedDim_facet <- function(type = "TSNE", facet_by){
  
  main_plot <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=round)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    labs(title = pca) +
    theme(legend.position = "none")
  
  facet_plot <- main_plot + facet_wrap(~round)
    
  return(main_plot + facet_plot)  
  
}
