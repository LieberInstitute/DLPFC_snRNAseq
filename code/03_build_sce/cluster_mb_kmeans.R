## adapted from https://www.stephaniehicks.com/biocdemo/articles/Demo.html

library("SingleCellExperiment")
library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")
library("numform")
library("jaffelab")

source("utils.R")
## load data
load(here("processed-data", "03_build_sce","sce_harmony_Sample.Rdata"), verbose = TRUE)

## neron paper returned 19 clusters for DLPFC, try 5:50
set.seed(610)
k_list <- seq(5, 50)

message("Apply mbkmeans from 5:50 - ", Sys.time())
km_res <- lapply(k_list, function(k) {
  message("k=",k)
  mbkmeans(sce, clusters = k, 
           batch_size = 500,
           reduceMethod = "HARMONY",
           calc_wcss = TRUE)
})

names(km_res[[1]])

save(km_res, file = here("processed-data", "03_build_sce","km_res.Rdata"))
# load(here("processed-data", "03_build_sce","km_res.Rdata"),verbose = TRUE)

## Use fasthplus to refine k
# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
l = lapply(km_res, `[[` ,5)
length(l)

find_t <- function(L, proportion = 0.05) {
  initial_t <- floor(length(L) * proportion)
  # message("init. t = ", initial_t)
  smallest_cluster_size <- min(table(L))
  n_labels <- length(unique(L))
  message("smallest cluster: ", smallest_cluster_size, ", n lables: ", n_labels)
  ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

message("Find fasthplus for clusters - ", Sys.time())
fasthplus <- lapply(l[20:46], function(li) {
  message(Sys.time())
  initial_t <- find_t(L=li,proportion = 0.01)
  h <- hpb(D=reducedDims(sce)$HARMONY,
           L=li,
           t=initial_t,
           r=30)
})

fasthplus <- unlist(fasthplus)

wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))
r <- read.csv(here("processed-data","03_build_sce","mb_kmeans_metrics.csv"))
colnames(r)[4:49]

results <- data.frame (k=k_list, wcss = wcss, fasthplus = fasthplus)
write.csv(results,file = here("processed-data","03_build_sce","mb_kmeans_metrics.csv"))

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_wcss.pdf"))
plot(k_list, wcss, type = "b")
abline(v=29, lty=2, col="red")
dev.off()


pdf(here("plots","03_build_sce","cluster", "mb_kmeans_fasthplus.pdf"))
plot(k_list, fasthplus, type = "b")
abline(v=29, lty=2, col="red")
dev.off()

#### Define SCE clusters with chosen k ####
sce$kmeans<- as.factor(paste0("mbk", f_pad_zero(km_res[[which(k_list==29)]]$Clusters)))

(cluster_tab <- table(sce$kmeans))
# mbk01 mbk02 mbk03 mbk04 mbk05 mbk06 mbk07 mbk08 mbk09 mbk10 mbk11 mbk12 mbk13 mbk14 mbk15 mbk16 mbk17 mbk18 mbk19 mbk20 
# 4233   264  1420 28963    35     2  2791   733  1235  1146  1192    16  4659   174  1587  3542   466  1645  1518   127 
# mbk21 mbk22 mbk23 mbk24 mbk25 mbk26 mbk27 mbk28 mbk29 
# 536  2718  1427  1615  5900   113  1193  2235  6119

summary(as.numeric(table(sce$kmeans)))

table(sce$kmeans, sce$round)
table(sce$kmeans, sce$Sample)
table(sce$kmeans, sce$subject)

#### check doublet score for each prelim clust ####
clusIndexes = splitit(sce$kmeans)
prelimCluster.medianDoublet <- sapply(clusIndexes, function(ii){
  median(sce$doubletScore[ii])
}
)

summary(prelimCluster.medianDoublet)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1370  0.4148  0.6658  0.6797  0.8024  1.5806 


## Plotting 

# Plotting set up 
cluster_colors <- DeconvoBuddies::create_cell_colors(cell_types = levels(sce$kmeans), pallet = "gg")

my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","cluster")

## Plot clusters in TSNE
TSNE_clusters <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=kmeans)) +
  geom_point(size = 0.2, alpha = 0.3) +
  scale_color_manual(values = cluster_colors) +
  my_theme +
  coord_equal()

ggsave(TSNE_clusters + 
         guides(colour = guide_legend(override.aes = list(size=2))),
       filename = here(plot_dir, "clusters_mbkm-29.png"), width = 10)

ggsave(TSNE_clusters + 
         facet_wrap(~kmeans) + 
         theme(legend.position = "none"), 
       filename = here(plot_dir, "clusters_mbkm-29_facet.png"), width = 10, height = 10)


#### mb kmeans Vs. walk trap ####
load(here("processed-data", "03_build_sce", "clusters.Rdata"), verbose = TRUE)
cvk <- table(clusters, sce$kmeans)

library(pheatmap)

pdf(here(plot_dir, "cluster_pheat_test.pdf"))
pheatmap(cvk)
dev.off()


#### Marker Genes ####
# Just for logcounts
sce <- batchelor::multiBatchNorm(sce, batch=sce$Sample)

# load Mathy's markers
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/revision/markers.rda", verbose = TRUE)
# markers.mathys.custom
all(unlist(markers.mathys.custom) %in% rowData(sce)$gene_name)

rownames(sce) <- rowData(sce)$gene_name

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_29_mathys_markers.pdf"), height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  message(names(markers.mathys.custom)[[i]])
  print(
    plotExpressionCustom(sce = sce,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "kmeans")+
      scale_color_manual(values = cluster_colors)
  )
}
dev.off()

#### Tran Maynard Top Markers####
dlpfc_markers <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv")
dlpfc_markers_list <- as.list(dlpfc_markers[,grepl("_1vAll", colnames(dlpfc_markers))])

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_29_Tran_markers.pdf"), height=6, width=8)
for(i in 1:length(dlpfc_markers_list)){
  message(names(dlpfc_markers_list)[[i]])
  f <- dlpfc_markers_list[[i]]
  f_good <- f[f %in% rownames(sce)]
  if(f != f_good) message("Missing...",paste(f[!f %in% f_good], collapse = ", "))
  print(
    plotExpressionCustom(sce = sce,
                         features = f_good,
                         features_name = names(dlpfc_markers_list)[[i]],
                         anno_name = "kmeans")+
      scale_color_manual(values = cluster_colors)
  )
}
dev.off()

#### Mean Ratio Top Markers####

load("dcl01/lieber/ajaffe/lab/deconvolution_bsp2/data/marker_stats_pan.v2.Rdata", verbose = TRUE)

# dlpfc_markers <- read.csv("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/tables/revision/top40genesLists_DLPFC-n3_cellType_SN-LEVEL-tests_LAH2020.csv")
# dlpfc_markers_list <- as.list(dlpfc_markers[,grepl("_1vAll", colnames(dlpfc_markers))])

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_29_Tran_markers.pdf"), height=6, width=8)
for(i in 1:length(dlpfc_markers_list)){
  message(names(dlpfc_markers_list)[[i]])
  f <- dlpfc_markers_list[[i]]
  f_good <- f[f %in% rownames(sce)]
  if(f != f_good) message("Missing...",paste(f[!f %in% f_good], collapse = ", "))
  print(
    plotExpressionCustom(sce = sce,
                         features = f_good,
                         features_name = names(dlpfc_markers_list)[[i]],
                         anno_name = "kmeans")+
      scale_color_manual(values = cluster_colors)
  )
}
dev.off()

#### Test Plots ####

test_plot <- ggcells(sce, mapping=aes(x=kmeans, y=SLC17A7, fill = kmeans)) + 
  geom_violin(scale = "width") +
  theme(legend.position = "None") +
  theme_bw()

ggsave(test_plot, filename = here("plots","03_build_sce","cluster", "test_ggcells_exprs.png"), width = 15)


## note small clusters
(small_clusters <- cluster_tab[cluster_tab < 20])
# mbk05 mbk18 mbk24 mbk29 
# 11     3     2     7

#### Annotate with marker cell types ####
anno <- read.csv(here("processed-data", "03_build_sce", "DLPFC_k29_anno.csv"))
table(anno$broad)

sce$cellType.broad <- factor(anno$broad[match(sce$kmeans, anno$cluster)])
table(sce$cellType.broad)
# Astro Excit Inhib Micro Nural Oligo   OPC small 
# 3557 18035 11380  2242  1330 39257  1791    12 

## Precentage cell composition
100*round(table(sce$cellType.broad)/ncol(sce),3)
# Astro Excit Inhib Micro Nural Oligo   OPC small 
# 4.6  23.2  14.7   2.9   1.7  50.6   2.3   0.0

load("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/TREG_paper/processed-data/00_data_prep/cell_colors.Rdata", verbose = TRUE)


pdf(here("plots","03_build_sce","cluster", "mb_broad_mathys_markers.pdf"), height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  message(names(markers.mathys.custom)[[i]])
  print(
    plotExpressionCustom(sce = sce,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "cellType.broad")+
      scale_color_manual(values = cell_colors)
  )
}
dev.off()



# sgejobs::job_single('cluster_mb_kmeans', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript cluster_mb_kmeans.R")
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
