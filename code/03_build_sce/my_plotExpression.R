
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
my_plotExpression <- function(sce, genes, assay = "logcounts", cat = "cellType", fill_colors = NULL){
  
  cat_df <- as.data.frame(colData(sce))[,cat, drop = FALSE]
  expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes,])) 
  
  cat <- cat_df[expression_long$Var2,]
  expression_long <- cbind(expression_long, cat)

  expression_violin <- ggplot(data = expression_long, aes(x = cat, y = value, fill = cat)) +
    geom_violin() +
    facet_wrap(~Var1)+
    labs(y = paste0("Expression (", assay,")")) +
    theme_bw() +
    theme(legend.position = "None",axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90,hjust=1))
  
  if(!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)

  # expression_violin
  return(expression_violin)
  
}

