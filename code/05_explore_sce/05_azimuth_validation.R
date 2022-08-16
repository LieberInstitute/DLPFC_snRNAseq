######### Cell Annotation with Azimuth#####
# devtools::install_github("satijalab/seurat-data") # install SeuratData package
# devtools::install_github("satijalab/azimuth", ref = "release/0.4.5") # install Azimuth package

library("SingleCellExperiment")
library("Seurat")
library("Azimuth")
library("SeuratData")
library("patchwork")
library("here")
library("pheatmap")


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


## Cell Annotation with Azimuth for each sample ####

## Use Azimuth Refrence ##
## new data
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

# create a Seurat object for each sample from SingleCellExperiment object
# query <- lapply(singlet.sce, function(sce)  CreateSeuratObject(counts = counts(sce), 
#                                                                meta.data = data.frame(colData(sce)),
#                                                                project = "SD"))

query <- CreateSeuratObject(counts = as.matrix(counts(sce)), 
                            meta.data = data.frame(colData(sce)),
                            project = "DLPFC")

rm(sce)

## Save 
# SeuratDisk::SaveH5Seurat(query, filename = here("processed-data","05_explore_sce","05_azimuth_validation","sce_DLPFC.h5Seurat"))

# query <- RunAzimuth(query, reference = ref.dir) ## Cell annotation with Azimuth

query <- RunAzimuth(query, reference = "humancortexref") ## Cell annotation with Azimuth

table(query$predicted.subclass)
# Astro       Endo    L2/3 IT      L5 ET      L5 IT    L5/6 NP      L6 CT      L6 IT L6 IT Car3        L6b 
# 6744       9888      21440        109       6775        345       1158       1155        332        875 
# Lamp5  Micro-PVM      Oligo        OPC      Pvalb       Sncg        Sst  Sst Chodl        Vip       VLMC 
# 802       2507      11384       1836       1513        243       4689         87       4178       1544 

ct_annos <- table(query$cellType_hc, query$predicted.subclass)


ct_annos_prop <- sweep(ct_annos,1,rowSums(ct_annos),`/`)
rowSums(ct_annos_prop)

ct_max <- apply(ct_annos_prop, 1, max)
table(ct_max  > 0.7)

sum(ct_annos_prop > .90)

png(here(plot_dir, "azimuth_v_hc.png"),height = 800, width = 800)
pheatmap(ct_annos_prop
         # ,
         # annotation_col= hc_anno,
         # annotation_row = km_anno,
         # annotation_colors = list(broad = temp_pallet)
)
dev.off()

library(ggplot2)
p1 <- DimPlot(query, group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
p2 <- DimPlot(query, group.by = "cellType_hc", label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
p3 <- DimPlot(query, group.by = "cellType_broad_hc", label = TRUE, label.size = 3, repel = TRUE) +
  NoLegend()
p4 <- DimPlot(query, group.by = "Sample")

ggsave((p1 + p2)/(p3 + p4), filename = here(plot_dir, "Azimuth_UMAP.png"), width = 13, height = 10)


# query <- lapply(query, function(so) RunAzimuth(so, reference = "reference/")) ## Cell annotation with Azimuth

# In predicted.subclass there are labels obtained from the annotation

## Integration####
query <- NormalizeData(query, verbose = TRUE)
# query <- lapply(query, NormalizeData, verbose = TRUE) # normalization for each sample

query <- FindVariableFeatures(query, nfeatures = 2e3, selection.method = "vst", verbose = TRUE)

# query <- lapply(query, FindVariableFeatures, nfeatures = 2e3,
#                 selection.method = "vst", verbose = TRUE) # find most variable features for each sample

query <- ScaleData(query, verbose = TRUE)
# query <- lapply(query, ScaleData, verbose = TRUE) # scale the data, for each sample

# find anchors & integrate
as <- FindIntegrationAnchors(query, verbose = FALSE)
seurat.obj.integrated <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE) # Seurat object

DefaultAssay(seurat.obj.integrated) <- "integrated"
# scale integrated data
seurat.obj.integrated <- ScaleData(seurat.obj.integrated, verbose = FALSE) # rescaled 

seurat.obj.integrated <- RunPCA(seurat.obj.integrated, assay = "integrated", npcs = 30, verbose = FALSE,
                                reduction.name = "PCA")
seurat.obj.integrated <- RunTSNE(seurat.obj.integrated, assay = "integrated", reduction = "PCA", dims = seq_len(20),
                                 seed.use = 1, do.fast = TRUE, verbose = FALSE,
                                 reduction.name = "TSNE")
seurat.obj.integrated <- RunUMAP(seurat.obj.integrated, assay = "integrated", reduction = "PCA", dims = seq_len(20),
                                 seed.use = 1, do.fast = TRUE, verbose = FALSE,
                                 reduction.name = "UMAP")


saveRDS(object = seurat.obj.integrate, file = file.path(ref.dir, "seurat.obj.integrate.Rds"))

# sgejobs::job_single('05_azimuth_validation', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript 05_azimuth_validation.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

