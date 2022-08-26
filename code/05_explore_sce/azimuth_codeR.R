######### Cell Annotation with Azimuth#####
devtools::install_github("satijalab/seurat-data") # install SeuratData package
devtools::install_github("satijalab/azimuth", ref = "release/0.4.5") # install Azimuth package

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(here)

# transform reference dataset to a Seurat object compatible with Azimuth####
## What refrence should we use?
# reference = "humancortexref"
# reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")
# Can we use Azimuth style refrence? or do we need to create a new one?

reference.so <- CreateSeuratObject(counts = counts(reference), meta.data = data.frame(colData(reference))) # where reference is a SingleCellExperiment object

## Normalize data
reference.so <- SCTransform(reference.so,
    assay = "RNA", new.assay.name = "SCT", variable.features.n = 2000,
    verbose = TRUE, conserve.memory = TRUE
) # conserve.memory=TRUE for large dataset (no R crash)

reference.so <- RunPCA(reference.so,
    assay = "SCT", npcs = 50, verbose = FALSE,
    reduction.name = "PCA", return.model = TRUE
) # npcs must be 50
reference.so <- RunUMAP(reference.so,
    assay = "SCT", reduction = "PCA", dims = seq_len(50),
    seed.use = 1, verbose = FALSE, reduction.name = "umap", return.model = TRUE
) # reduction.name = "umap" to be able to create the object of interest

reference.so$subclass <- as.factor(reference.so$subclass)
Idents(object = reference.so) <- "subclass"


reference.azimuth <- AzimuthReference(reference.so,
    refUMAP = "umap", refDR = "PCA", refAssay = "SCT", dims = 1:50,
    metadata = c("subclass"), verbose = TRUE
) # object compatible with Azimuth

# save reference
ref.dir <- here("processed-data", "05_explore_sce", "05_azimuth_validataion", "reference")

SaveAnnoyIndex(object = reference.azimuth[["refdr.annoy.neighbors"]], file = file.path(ref.dir, "idx.annoy"))
saveRDS(object = reference.azimuth, file = file.path(ref.dir, "ref.Rds"))

## Cell Annotation with Azimuth for each sample ####
# create a Seurat object for each sample from SingleCellExperiment object

## should we also loop by Sample?
query <- lapply(singlet.sce, function(sce) {
    CreateSeuratObject(
        counts = counts(sce),
        meta.data = data.frame(colData(sce)),
        project = "SD"
    )
})

query <- lapply(query, function(so) RunAzimuth(so, reference = ref.dir)) ## Cell annotation with Azimuth
# In predicted.subclass there are labels obtained from the annotation

## What happens past this point? Is it relevant for cell type annotation?

## Integration####
## why not SCTransform?
query <- lapply(query, NormalizeData, verbose = TRUE) # normalization for each sample
query <- lapply(query, FindVariableFeatures,
    nfeatures = 2e3,
    selection.method = "vst", verbose = TRUE
) # find most variable features for each sample
query <- lapply(query, ScaleData, verbose = TRUE) # scale the data, for each sample

# find anchors & integrate
as <- FindIntegrationAnchors(query, verbose = FALSE)
seurat.obj.integrated <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE) # Seurat object
DefaultAssay(seurat.obj.integrated) <- "integrated"

# scale integrated data
seurat.obj.integrated <- ScaleData(seurat.obj.integrated, verbose = FALSE) # rescaled
seurat.obj.integrated <- RunPCA(seurat.obj.integrated,
    assay = "integrated", npcs = 30, verbose = FALSE,
    reduction.name = "PCA"
)
seurat.obj.integrated <- RunTSNE(seurat.obj.integrated,
    assay = "integrated", reduction = "PCA", dims = seq_len(20),
    seed.use = 1, do.fast = TRUE, verbose = FALSE,
    reduction.name = "TSNE"
)
seurat.obj.integrated <- RunUMAP(seurat.obj.integrated,
    assay = "integrated", reduction = "PCA", dims = seq_len(20),
    seed.use = 1, do.fast = TRUE, verbose = FALSE,
    reduction.name = "UMAP"
)
