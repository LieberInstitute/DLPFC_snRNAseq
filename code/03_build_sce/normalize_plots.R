library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")

# Plotting set up 
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","normalize2")

if(!dir.exists(plot_dir)) dir.create(plot_dir)

## List normalized files
(norm_fn <- list.files(path = here("processed-data", "03_build_sce"), pattern = "^sce",
                      full.names = TRUE))

names(norm_fn) <- gsub(".Rdata", "", gsub("sce_","",basename(norm_fn)))

load(norm_fn[[1]], verbose = TRUE)

walk2(norm_fn, names(norm_fn), function(sce_fn, name){
  
  annotation <- gsub("_", " + ", name)
  message(annotation)
  ## Load 
  load(sce_fn, verbose = TRUE)
  all(c("TSNE", "UMAP") %in% reducedDimNames(sce))
  
  #### Plot TSNE ####
  
  ## TSNE_method-correction_color-colorvar.png
  ## ex. TSNE_MNN-round_color-subject.png
  
  TSNE_plot_name = paste0("TSNE_", gsub("PCA_","",pca))
  
  ## ggplots 
  tsne_round <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=round)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    labs(title = pca)
  
  ggsave(tsne_round, filename = here(plot_dir, paste0(plot_name,  ".png")), width = 10)
  ggsave(tsne_round +
           theme(legend.position = "none") + 
           tsne_round + facet_wrap(~round) +
           theme(legend.position = "none"), 
         filename = here(plot_dir,paste0(plot_name, "_round.png")), width = 20, height = 10)
  
  
  # ggsave(tsne_round + facet_wrap(~Sample), filename = here(plot_dir,"MNN_TSNE_roundXsample.png"), width = 10)
  
  # ggsave(tsne_round + facet_grid(round~subject), filename = here(plot_dir,"sce_TSNE_roundXsubject.png"), width = 10)
  
  
  tsne_sample <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=Sample)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(tsne_sample + tsne_sample + facet_wrap(~Sample), 
         filename = here(plot_dir,paste0(plot_name, "_sample.png")),
         width = 20, height = 10)
  
  tsne_subject <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=subject)) +
    geom_point(size = 0.2, alpha = 0.3) +
    coord_equal() +
    theme_bw() +
    theme(legend.position = "none")
  
  ggsave(tsne_subject + tsne_subject + facet_wrap(~subject), 
         filename = here(plot_dir, paste0(plot_name,"_subject.png")),
         width = 20, height = 10)
  
  
  # tsne_discard <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=discard_auto)) +
  #   geom_point(size = 0.2, alpha = 0.3) +
  #   coord_equal() +
  #   theme_bw() 
  # 
  # ggsave(tsne_discard, filename = here(plot_dir,"MNN_TSNE_discard_auto.png"))
  # 
  # tsne_doubletScore <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=doubletScore)) +
  #   geom_point(size = 0.2, alpha = 0.3) +
  #   coord_equal() +
  #   theme_bw() +
  #   viridis::scale_color_viridis(trans = "log")
  # 
  # ggsave(tsne_doubletScore, filename = here(plot_dir,"MNN_TSNE_doubletScore.png"))
  
})

## Save data
save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('normalize', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
