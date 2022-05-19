library("SingleCellExperiment")
library("scran")
library("scater")
# # library("jaffelab")
# library("tidyverse")
library("batchelor")
library("here")
library("sessioninfo")

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
uncorrected <- runPCA(sce, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())

uncorrected <- runTSNE(uncorrected, dimred="PCA")

save(uncorrected, file = here("processed-data", "03_build_sce","uncorrected_TSNE.Rdata"))

plotTSNE(uncorrected, colour_by="batch")

#### Preform Batch Correction with MNN ####

## should we use "Sample" (19), "Round" (5), or "subject" (?) for correction?

table(sce$round)
# round0 round1 round2 round3 round4 round5 
# 3426  10797  16671  20606  16509  16747 

# Run `fastMNN` (internally uses `multiBatchPCA`), taking top 100 (instead of default 50 PCs)
set.seed(519)
mnn.hold <-  fastMNN(sce, batch=sce$round,
                     subset.row=chosen.hvgs, d=100,
                     correct.all=TRUE, get.variance=TRUE,
                     BSPARAM=BiocSingular::IrlbaParam())

# Add them to the SCE, as well as the metadata
reducedDim(sce, "PCA_corrected") <- reducedDim(mnn.hold, "corrected") # 100 components
metadata(sce) <- metadata(mnn.hold)

## Should we find optimal PC space? - use all 100 for now
## t-SNE
set.seed(109)
sce <- runTSNE(sce, dimred="PCA_corrected")


## UMAP
set.seed(109)
sce <- runUMAP(sce, dimred="PCA_corrected")


# How do these look?
pdf(file = here("plots","03_build_sce","MNN_TSNE.pdf"))
plotReducedDim(sce, dimred="TSNE", colour_by="round")
plotReducedDim(sce, dimred="TSNE", colour_by="subject")
plotReducedDim(sce, dimred="TSNE", colour_by="Sample")
dev.off()

pdf(file = here("plots","03_build_sce","MNN_UMAP.pdf"))
plotReducedDim(sce, dimred="UMAP", colour_by="round")
plotReducedDim(sce, dimred="UMAP", colour_by="subject")
plotReducedDim(sce, dimred="UMAP", colour_by="Sample")
dev.off()

plotReducedDim(sce, dimred="UMAP", colour_by="sampleID")

## Save data
save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))

# sgejobs::job_single('normalize_and_doublet_score', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize_and_doublet_score.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
