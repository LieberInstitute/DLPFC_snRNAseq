library("SingleCellExperiment")
library("spatialLIBD")
library("tidyverse")
library("scMerge")
library("here")
library("sessioninfo")

## Plot setup
plot_dir <- here("plots", "05_explore_sce", "09_SCZ_interaction_expression")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# data_dir <- here("processed-data", "05_explore_sce", "09_SCZ_interaction_expression")
# if (!dir.exists(data_dir)) dir.create(data_dir)

## Load data
load(here("processed-data", "05_explore_sce", "08_pseudobulk_cellTypes", "sce_pseudo-cellType_layer.Rdata"), verbose = TRUE)
## Check out colData
table(sce_pseudo$cellType_layer)
# Astro    EndoMural        Micro        Oligo          OPC   Excit_L2/3     Excit_L3 Excit_L3/4/5     Excit_L4 
# 18           17           18           19           19            1           19           17           16 
# Excit_L5   Excit_L5/6     Excit_L6        Inhib 
# 16           16           16           19


## Genes to explore
SCZ_interactions <- read.csv("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/rdata/spe/14_spatial_registration_PEC/04_PEC_check_expression/SCZ_top_interactions.csv")
genes_of_interest <- unique(unlist(SCZ_interactions))

genes_of_interest[!genes_of_interest %in% rownames(sce_pseudo)]
genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(sce_pseudo)]
length(genes_of_interest)

#### PLOT Expression ####
cat_df <- as.data.frame(colData(sce_pseudo))[, c("SAMPLE_ID", "cellType_layer"), drop = FALSE]

expression_long <- reshape2::melt(as.matrix(logcounts(sce_pseudo)[genes_of_interest, , drop = FALSE]))
cat <- cat_df[expression_long$Var2, ]
colnames(cat) <- c("sample_id", "cat")
expression_long <- cbind(expression_long, cat)

pdf(here(plot_dir, "spatialDLPFC_SCZ_interactions.pdf"))
for(i in 1:nrow(SCZ_interactions)){
  
  interaction <- SCZ_interactions[i,]
  interaction_title <- paste(interaction[[1]], "->", interaction[[2]])
  message(i, " ", interaction_title)
  
  expression_temp <- expression_long |>
    filter(Var1 %in% interaction) |>
    mutate(Var1 = factor(Var1, levels = interaction))
  
  expression_violin <- ggplot(data = expression_temp, aes(x = cat, y = value, fill = cat)) +
    geom_violin(scale = "width", alpha = 0.5) +
    geom_jitter(aes(color = cat), size = .3) +
    facet_wrap(~Var1, ncol = 1, scales = "free_y") +
    labs(
      title = interaction_title,
      y = "Expression (logcounts)"
    ) +
    theme_bw() +
    theme(
      legend.position = "None",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text.x = element_text(face = "italic"),
      text = element_text(size = 15)
    ) +
    stat_summary(
      fun = median,
      # fun.min = median,
      # fun.max = median,
      geom = "crossbar",
      width = 0.3
    ) +
    scale_color_manual(values = metadata(sce_pseudo)$cell_type_colors_layer) +
    scale_fill_manual(values = metadata(sce_pseudo)$cell_type_colors_layer)
  print(expression_violin)
}
dev.off()


