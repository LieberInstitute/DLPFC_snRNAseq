######### Cell Annotation with Azimuth#####
# devtools::install_github("satijalab/seurat-data") # install SeuratData package
# devtools::install_github("satijalab/azimuth", ref = "release/0.4.5") # install Azimuth package

library("SingleCellExperiment")
library("Seurat")
library("Azimuth")
library("SeuratData")
# library("patchwork")
library("here")
# library("pheatmap")


#### Plotting ####
plot_dir = here("plots","05_explore_sce","05_azimuth_validation")
if(!dir.exists(plot_dir)) dir.create(plot_dir)


## Use Neuron DLPFC as refrence?
# load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/regionSpecific_DLPFC-n3_cleaned-combined_SCE_LAH2021.rda", verbose = TRUE)

# transform reference dataset to a Seurat object compatible with Azimuth####

## Only use top 2k deviant genes
# hdgs.hb <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing=T)][1:2000]

reference.so <- CreateSeuratObject(counts = as.matrix(counts(sce.dlpfc)), meta.data = data.frame(colData(sce.dlpfc))) # where reference is a SingleCellExperiment object
# reference.so <- CreateSeuratObject(counts = as.matrix(counts(sce[hdgs.hb,])), meta.data = data.frame(colData(sce))) # where reference is a SingleCellExperiment object
# reference.so <- CreateSeuratObject(counts = counts(reference), meta.data = data.frame(colData(reference))) # where reference is a SingleCellExperiment object

reference.so <- SCTransform(reference.so,assay = "RNA", new.assay.name = "SCT", variable.features.n = 2000,
                            verbose = TRUE,conserve.memory=TRUE) # conserve.memory=TRUE for large dataset (no R crash)

reference.so <- RunPCA(reference.so, assay = "SCT", npcs = 50, verbose = FALSE,
                       reduction.name = "PCA",return.model=TRUE) # npcs must be 50

reference.so <- RunUMAP(reference.so, assay = "SCT", reduction = "PCA", dims = seq_len(50),
                        seed.use = 1, verbose = FALSE, reduction.name = "umap",return.model=TRUE) # reduction.name = "umap" to be able to create the object of interest

reference.so$subclass <- as.factor(reference.so$cellType) 
Idents(object = reference.so) <- "subclass"

reference.azimuth <- AzimuthReference(reference.so, refUMAP = "umap",refDR = "PCA",refAssay = "SCT",dims = 1:50,
                                      metadata = c("subclass"), verbose = TRUE) # object compatible with Azimuth

# save reference
ref.dir <- here("processed-data","05_explore_sce","05_azimuth_validation","reference")
SaveAnnoyIndex(object = reference.azimuth[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = reference.azimuth, file = file.path(ref.dir, "ref.Rds"))


saveRDS(object = seurat.obj.integrate, file = file.path(ref.dir, "seurat.obj.integrate.Rds"))

# sgejobs::job_single('05_azimuth_validation', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript 05_azimuth_validation.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

