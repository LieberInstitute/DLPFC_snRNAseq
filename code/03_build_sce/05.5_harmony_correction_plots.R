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
my_theme <- theme_classic() +
    theme(text = element_text(size = 15))

plot_dir <- here("plots", "03_build_sce", "05_harmony_correction")

if (!dir.exists(plot_dir)) dir.create(plot_dir)

## List normalized files
sce_fn <- here("processed-data", "03_build_sce", c("sce_uncorrected_glm.Rdata", "sce_harmony_Sample.Rdata"))
all(file.exists(sce_fn))

# (sce_fn <- list.files(path = here("processed-data", "03_build_sce"), pattern = "^sce",
#                       full.names = TRUE))

sce_fn <- here("processed-data", "03_build_sce",c("sce_uncorrected_glm.Rdata", "sce_harmony_Sample.Rdata"))

(names(sce_fn) <- gsub(".Rdata", "", gsub("sce_", "", basename(sce_fn))))


walk2(sce_fn, names(sce_fn), function(sce_fn, name) {
  annotation <- gsub("_", " + ", name)
  message(annotation)
  ## Load
  load(sce_fn, verbose = TRUE)
  
  if (grepl("uncorrected", name)) sce <- sce_uncorrected
  all(c("TSNE", "UMAP") %in% reducedDimNames(sce))
  
  for (t in c("TSNE", "UMAP")) {
    
    ## Catagorical plots
    for (cat in c("round", "subject", "Sample")) {
      fn <- paste0(t, "_", name, "-color_", cat)
      
      message("plotting... ", fn)
      cat_plot <- plot_reducedDim_facet(sce, type = t, facet_by = cat)
      ggsave(cat_plot, filename = here(plot_dir, paste(fn,".png")), width = 13)
      ggsave(cat_plot, filename = here(plot_dir, paste(fn,".pdf")), width = 13)
      
    }
  }
})

# sgejobs::job_single('05.5_harmony_correction_plots', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 05.5_harmony_correction_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
