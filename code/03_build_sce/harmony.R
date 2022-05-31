library("SingleCellExperiment")
# library("scran")
library("scater")
library("batchelor")
library("harmony")
library("here")
library("sessioninfo")
library("patchwork")
library("ggplot2")

# Plotting themes
# my_theme <- theme_bw() +
#   theme(text = element_text(size=15))

## Load empty-free sce data
load(here("processed-data", "sce", "sce_no_empty_droplets.Rdata"), verbose = TRUE)

load(here("processed-data", "03_build_sce","uncorrected_TSNE.Rdata"), verbose = TRUE)

sce <- sce[,!sce$discard_auto]
dim(sce)
# [1] 36601 77604

sce <- multiBatchNorm(sce, batch=sce$round)

sce34 <- sce[,sce$round %in% c("round3", "round4")]

# V <- cell_lines$scaled_pcs
# V <- assays(sce34)$counts
V  <- reducedDims(uncorrected)$PCA
class(V)

# meta_data <- cell_lines$meta_data
meta_data <- as.data.frame(colData(uncorrected))

harmony_embeddings <- harmony::HarmonyMatrix(
  as.matrix(V), meta_data, 'round', do_pca = FALSE, verbose=FALSE
)

sce2 <- RunHarmony(sce, 'round', lambda = .1, verbose = TRUE)

V2 <- as.data.frame(cbind(V,meta_data[,"round",drop = FALSE]))
harmony_embeddings2 <- as.data.frame(cbind(harmony_embeddings,meta_data[,"round",drop = FALSE]))

before <- ggplot(data = V2, aes(x = PC1, y = PC2, color = round)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw()

after_harmony <- ggplot(harmony_embeddings2, aes(x = PC1, y = PC2, color = round)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw()

ggsave(before + after_harmony, filename = here("plots","03_build_sce","normalize", "harmony_test.png"), width = 14)

sce <- uncorrected
reducedDim(sce, "PCA_harmony") <- harmony_embeddings

set.seed(109)
sce <- runTSNE(sce, dimred = "PCA_harmony")
sce <- runUMAP(sce, dimred = "PCA_harmony")

#### TSNE Plots ####

plot_name = paste0("Harmony_TSNE_", "round")

tsne_round <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=round)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  theme_bw() +
  labs(title = pca)

ggsave(tsne_round, filename = here(plot_dir, paste0(plot_name,  ".png")), width = 10)

