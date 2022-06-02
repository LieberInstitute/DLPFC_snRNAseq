library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")
library("purrr")

source(here("code", "03_build_sce", "utils.R"))

# Plotting set up 
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","normalize2")

if(!dir.exists(plot_dir)) dir.create(plot_dir)

## List normalized files
(norm_fn <- list.files(path = here("processed-data", "03_build_sce"), pattern = "^sce",
                      full.names = TRUE))

names(norm_fn) <- gsub(".Rdata", "", gsub("sce_","",basename(norm_fn)))

# load(norm_fn[[1]], verbose = TRUE)
test_plot <- plot_reducedDim_qc(sce, title = "MNN + round")
ggsave(test_plot, filename = here(plot_dir,"TSNE_test.png"),
       width = 10, height = 10)
# 


walk2(norm_fn[3], names(norm_fn)[3], function(sce_fn, name){
  annotation <- gsub("_", " + ", name)
  message(annotation)
  ## Load 
  load(sce_fn, verbose = TRUE)
  
  if(name == "uncorrected") sce <- sce_uncorrected
  all(c("TSNE", "UMAP") %in% reducedDimNames(sce))
  
  for(t in c("TSNE","UMAP")){

    ## Catagorical plots
    for(cat in c("round","subject")){

      fn = paste0(t,"_",name,"-color_",cat,".png")
      message(fn)

      cat_plot <- plot_reducedDim_facet(sce, type = t, facet_by = cat, title = annotation)
      ggsave(cat_plot, filename = here(plot_dir,fn), width = 20, height = 10)

    }

  }
  
  # ## only TSNE for QC metrics plots
  # for(qc in c("sum","doubletScore")){
  #   fn = paste0("TSNE_",name,"-",qc,".png")
  #   message(fn)
  # 
  #   qc_plot <- plot_reducedDim_qc(sce, type = "TSNE", color_by = qc, title = annotation)
  #   ggsave(qc_plot, filename = here(plot_dir,fn), width = 10, height = 10)
  # }
  
})

# sgejobs::job_single('normalize', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
