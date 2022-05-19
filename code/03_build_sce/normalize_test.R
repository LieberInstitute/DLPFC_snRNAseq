library("SingleCellExperiment")
library("scran")
library("scater")
library("batchelor")
library("here")
library("sessioninfo")

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
pdf(file = here("plots","03_build_sce","MNN_TSNE.pdf"))
plotReducedDim(uncorrected, dimred="TSNE", colour_by="round")
plotReducedDim(uncorrected, dimred="TSNE", colour_by="subject")
plotReducedDim(uncorrected, dimred="TSNE", colour_by="Sample")
dev.off()

plotReducedDim(sce, dimred="UMAP", colour_by="sampleID")


# sgejobs::job_single('normalize_test', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript normalize_test.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
