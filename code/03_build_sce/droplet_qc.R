library("SingleCellExperiment")
library("DropletUtils")
library("tidyverse")
library("patchwork")
library("here")
library("sessioninfo")

## plot set up
my_theme <- theme_bw() +
  theme(text = element_text(size=15))


## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)
length(table(sce$Sample))

droplet_score_fn <- list.files(here("processed-data", "03_build_sce","droplet_scores"),
                               full.names = TRUE)

names(droplet_score_fn)  <- gsub("droplet_scores_|.Rdata","",basename(droplet_score_fn))

e.out <- lapply(droplet_score_fn, function(x) get(load(x)))
map_int(e.out, nrow)

## Check empty droplet results
map(e.out, ~addmargins(table(Signif = .x$FDR <= 0.001, Limited = .x$Limited, useNA = "ifany")))

pdf(here("plots","03_build_sce","droplet_qc_plots.pdf"))
walk2(names(e.out), e.out, function(sample, e){
  
  sce_temp <- sce[,sce$Sample == sample]
  br.out <- barcodeRanks(sce_temp)
  
  n_cell_anno <- paste("Non-empty:", sum(e$FDR < 0.01, na.rm = TRUE))
  message(sample, n_cell_anno)
  
  droplet_elbow_plot <- as.data.frame(br.out) %>%
    add_column(FDR = e$FDR) %>%
    ggplot(aes(x = rank, y = total, color = FDR < 0.01)) +
    geom_point(alpha = 0.5, size = 1) + 
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(x = "Barcode Rank", y = "Total UMIs",
         title = paste("Sample", sample), subtitle = n_cell_anno) +
    my_theme +
    theme(legend.position = "None")
  
  droplet_scatter_plot <- as.data.frame(e.out$Br2720_mid) %>%
    ggplot(aes(x = Total, y = -LogProb, color = FDR < 0.01)) +
    geom_point(alpha = 0.5, size = 1) + 
    labs(x = "Total UMIs", y = "-Log Probability") +
    my_theme+
    theme(legend.position = "bottom")
  print(droplet_elbow_plot/droplet_scatter_plot)
  ggsave(droplet_elbow_plot/droplet_scatter_plot, filename = here("plots","03_build_sce", "droplet_qc_png",paste0("droplet_qc_",sample,".png")))
  
})
dev.off()


## Eliminate empty droplets
sce <- sce[, which(e.out$FDR <= 0.001)]

## Compute QC metrics
sce <- addPerCellQC(
  sce,
  subsets = list(Mito = which(seqnames(sce) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

## Save for later
save(sce, file = here::here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

