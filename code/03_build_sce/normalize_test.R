library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")

## Load empty-free sce data
load(here("processed-data", "sce", "sce_no_empty_droplets.Rdata"), verbose = TRUE)

#### Rescale ####
# Use `multiBatchNorm()` to compute log-normalized counts
sce <- multiBatchNorm(sce, batch=sce$round)

# Find HVGs
geneVar <- modelGeneVar(sce)
chosen.hvgs <- geneVar$bio > 0
sum(chosen.hvgs)
# [1] 13254

#### How does the un-normalized data look? ####
uncorrected <- runPCA(sce, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())

uncorrected <- runTSNE(uncorrected, dimred="PCA")

save(uncorrected, file = here("processed-data", "03_build_sce","uncorrected_TSNE.Rdata"))

plotTSNE(uncorrected, colour_by="batch")

# How do these look?
pdf(file = here("plots","03_build_sce","normalize","uncorrected_TSNE.pdf"))
plotTSNE(uncorrected, colour_by="round")
plotTSNE(uncorrected, colour_by="subject")
plotTSNE(uncorrected, colour_by="Sample")
dev.off()

tsne_test <- ggcells(uncorrected, mapping=aes(x=TSNE.1, y=TSNE.2, colour=round)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw()

ggsave(tsne_test, filename = here("plots","03_build_sce","normalize","uncorrected_TSNE.png"), width = 10)
ggsave(tsne_test +
         theme(legend.position = "none") + 
         tsne_test + facet_wrap(~round) +
         theme(legend.position = "none"), 
       filename = here("plots","03_build_sce","normalize","uncorrected_TSNE_round.png"), width = 20, height = 10)


ggsave(tsne_test + facet_wrap(~Sample), filename = here("plots","03_build_sce","normalize","uncorrected_TSNE_roundXsample.png"), width = 10)

ggsave(tsne_test + facet_grid(round~subject), filename = here("plots","03_build_sce","normalize","uncorrected_TSNE_roundXsubject.png"), width = 10)


tsne_sample <- ggcells(uncorrected, mapping=aes(x=TSNE.1, y=TSNE.2, colour=Sample)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

ggsave(tsne_sample + tsne_sample + facet_wrap(~Sample), 
       filename = here("plots","03_build_sce","normalize","uncorrected_TSNE_sample.png"),
       width = 20, height = 10)

tsne_subject <- ggcells(uncorrected, mapping=aes(x=TSNE.1, y=TSNE.2, colour=subject)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

ggsave(tsne_subject + tsne_subject + facet_wrap(~subject), 
       filename = here("plots","03_build_sce","normalize","uncorrected_TSNE_subject.png"),
       width = 20, height = 10)

# sgejobs::job_single('normalize_test', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize_test.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
