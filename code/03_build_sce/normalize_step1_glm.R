library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("HDF5Array")
library("here")
library("sessioninfo")

## Load HD5F
sce <- loadHDF5SummarizedExperiment(dir=here("processed-data", "03_build_sce", "hdf5_sce"), prefix="")
dim(sce)

#### DEviance Feat Selection ####
set.seed(606)

message("running Deviance Feat. Selection - ", Sys.time())
sce <- devianceFeatureSelection(sce,
                                assay="counts", fam="binomial", sorted=F,
                                batch=as.factor(sce$round))

# This temp file just used for getting batch-corrected components (drops a variety of entries)

pdf(here("plots","03_build_sce","normalize1", "binomial_deviance.pdf"))
plot(sort(rowData(sce)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")
dev.off()

message("running nullResiduals - ", Sys.time())
sce <- nullResiduals(sce, assay="counts", fam="binomial",  # default params
                     type="deviance")


hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing=T)][1:2000]

message("running PCA - ", Sys.time())
sce_uncorrected <- runPCA(sce, exprs_values="binomial_deviance_residuals",
                      subset_row=hdgs.hb, ncomponents=100,
                      name="GLMPCA_approx",
                      BSPARAM=BiocSingular::IrlbaParam())

## Why multi match as well?
# sce <- multiBatchNorm(sce, batch=sce$sample_short)

#### Reduce Dims####

message("running TSNE - ", Sys.time())
sce_uncorrected <- runTSNE(sce_uncorrected, dimred="GLMPCA_approx")

message("running UMAP - ", Sys.time())
sce_uncorrected <- runUMAP(sce_uncorrected, dimred="GLMPCA_approx")

message("Saving Data - ", Sys.time())
save(sce_uncorrected, file = here("processed-data", "03_build_sce","sce_uncorrected_glm.Rdata"))

# sgejobs::job_single('normalize_step1_glm', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript normalize_step1_glm.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
