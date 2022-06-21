
#' Plot SCE Expression of multiple genes
#'
#' @param sce A SingleCellExperiment object containing expression values
#' @param genes list of genes to plot
#' @param assay name of assay to plot
#' @param cat name of category 
#' @param fill_colors optional color pallet for cat
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' my_plotExpression(example_sce, c("Gene_0001", "Gene_0002","Gene_0003", "Gene_0004"), cat = "Mutation_Status")
#' my_plotExpression(example_sce, 
#'                  c("Gene_0001", "Gene_0002","Gene_0003", "Gene_0004"), 
#'                  cat = "Mutation_Status",
#'                  fill_colors = c(negative = "green", positive = "pink"))
my_plotExpression <- function(sce, genes, assay = "logcounts", cat = "cellType", fill_colors = NULL, title = NULL){
  
  cat_df <- as.data.frame(colData(sce))[,cat, drop = FALSE]
  expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes,])) 
  
  cat <- cat_df[expression_long$Var2,]
  expression_long <- cbind(expression_long, cat)

  expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value, fill = cat)) +
    geom_violin(scale = "width") +
    facet_wrap(~Var1, ncol = 2)+
    labs(y = paste0("Expression (", assay,")"),
         title = title) +
    theme_bw() +
    theme(legend.position = "None",axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,hjust=1),
          strip.text.x = element_text(face = "italic")) +
    stat_summary(fun = median, 
                 # fun.min = median, 
                 # fun.max = median,
                 geom = "crossbar", 
                 width = 0.3)
  
  if(!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

  # expression_violin
  return(expression_violin)
  
}


my_plotMarkers <- function(sce, marker_list, assay = "logcounts", cat = "cellType", fill_colors = NULL, pdf_fn){
  message("plotting: ", pdf_fn)
  pdf(pdf_fn, height=6, width=8)
  
  for(i in 1:length(marker_list)){
    message(names(marker_list)[[i]])
    
    m <- marker_list[[i]]
    markers <- m[m %in% rownames(sce)]
    
    # if(m != markers) message("Missing...",paste(m[!m %in% markers], collapse = ", "))
    if(m != markers) message("Missing markers...")
    print(
      my_plotExpression(sce,
                        genes = markers, 
                        title = names(marker_list)[[i]],
                        cat = cat,
                        fill_colors = fill_colors)
    )
  }
  
  dev.off()
}

