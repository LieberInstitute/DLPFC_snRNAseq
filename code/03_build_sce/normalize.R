library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")
library("patchwork")

# Plotting themes
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

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
## move to test
# uncorrected <- runPCA(sce, subset_row=chosen.hvgs,
#                       BSPARAM=BiocSingular::RandomParam())
# 
# uncorrected <- runTSNE(uncorrected, dimred="PCA")
# 
# save(uncorrected, file = here("processed-data", "03_build_sce","uncorrected_TSNE.Rdata"))
# 
# plotTSNE(uncorrected, colour_by="batch")

#### Preform Batch Correction with MNN ####

## should we use "Sample" (19), "Round" (5), or "subject" (?) for correction?

table(sce$round)
# round0 round1 round2 round3 round4 round5 
# 3426  10797  16671  20606  16509  16747 

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(519)

## what about merge order? maybe use auto.merge=TRUE
message("running fast MNN", Sys.time())
mnn.hold <-  fastMNN(sce, batch=sce$round,
                     auto.merge=TRUE,
                     subset.row=chosen.hvgs, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())

message("Done MNN", Sys.time())

## Check out lost variance (looking for < 10%)
metadata(mnn.hold)$merge.info$lost.var
#           round0     round1      round2      round3      round4      round5
# [1,] 0.000000000 0.00000000 0.014235594 0.010002530 0.000000000 0.000000000
# [2,] 0.000000000 0.00000000 0.015129319 0.018261784 0.000000000 0.030320055
# [3,] 0.000000000 0.00000000 0.007791005 0.010618393 0.036963754 0.005388343
# [4,] 0.028918142 0.00000000 0.002539170 0.004024242 0.003849648 0.002530020
# [5,] 0.003774332 0.02697881 0.003331450 0.003361024 0.002730066 0.002931449

# Add them to the SCE, as well as the metadata
reducedDim(sce, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce) <- metadata(mnn.hold)

## Should we find optimal PC space? - use all 100 for now
## t-SNE
message("running TSNE + UMAP", Sys.time())

set.seed(109)
sce <- runTSNE(sce, dimred="PCA_corrected")
sce <- runUMAP(sce, dimred="PCA_corrected")

message("Done TSNE + UMAP", Sys.time())

# How do these look?
pdf(file = here("plots","03_build_sce","normalize","MNN_TSNE.pdf"))
plotReducedDim(sce, dimred="TSNE", colour_by="round")
plotReducedDim(sce, dimred="TSNE", colour_by="subject")
plotReducedDim(sce, dimred="TSNE", colour_by="Sample")
dev.off()

pdf(file = here("plots","03_build_sce","normalize","MNN_UMAP.pdf"))
plotReducedDim(sce, dimred="UMAP", colour_by="round")
plotReducedDim(sce, dimred="UMAP", colour_by="subject")
plotReducedDim(sce, dimred="UMAP", colour_by="Sample")
dev.off()

## Save data
save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

## ggplots 
tsne_test <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=round)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw()

ggsave(tsne_test, filename = here("plots","03_build_sce","normalize","MNN_TSNE.png"), width = 10)
ggsave(tsne_test +
         theme(legend.position = "none") + 
         tsne_test + facet_wrap(~round) +
         theme(legend.position = "none"), 
       filename = here("plots","03_build_sce","normalize","MNN_TSNE_round.png"), width = 20, height = 10)


ggsave(tsne_test + facet_wrap(~Sample), filename = here("plots","03_build_sce","normalize","MNN_TSNE_roundXsample.png"), width = 10)

ggsave(tsne_test + facet_grid(round~subject), filename = here("plots","03_build_sce","normalize","sce_TSNE_roundXsubject.png"), width = 10)


tsne_sample <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=Sample)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

ggsave(tsne_sample + tsne_sample + facet_wrap(~Sample), 
       filename = here("plots","03_build_sce","normalize","MNN_TSNE_sample.png"),
       width = 20, height = 10)

tsne_subject <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=subject)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")

ggsave(tsne_subject + tsne_subject + facet_wrap(~subject), 
       filename = here("plots","03_build_sce","normalize","MNN_TSNE_subject.png"),
       width = 20, height = 10)


tsne_discard <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=discard_auto)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() 

ggsave(tsne_discard, filename = here("plots","03_build_sce","normalize","MNN_TSNE_discard_auto.png"))

tsne_doubletScore <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=doubletScore)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  viridis::scale_color_viridis(trans = "log")

ggsave(tsne_doubletScore, filename = here("plots","03_build_sce","normalize","MNN_TSNE_doubletScore.png"))



# sgejobs::job_single('normalize', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
