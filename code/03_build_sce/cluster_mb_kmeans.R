## adapted from https://www.stephaniehicks.com/biocdemo/articles/Demo.html

library("SingleCellExperiment")
library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")

## load data
load(here("processed-data", "03_build_sce","sce_harmony_Sample.Rdata"), verbose = TRUE)

## neron paper returned 19 clusters for DLPFC, try 15:25
k_list <- seq(5, 50)

km_res <- lapply(k_list, function(k) {
  message("k=",k)
  mbkmeans(sce, clusters = k, 
           batch_size = 500,
           reduceMethod = "HARMONY",
           calc_wcss = TRUE)
})

wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_wcss.pdf"))
plot(k_list, wcss, type = "b")
abline(v=26, lty=2, col="red")
dev.off()

sce$kmeans <- paste0("mbk", km_res[[which(k_list==19)]]$Clusters)
table(sce$kmeans)

## Plotting 

# Plotting set up 
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","cluster")


ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=kmeans)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal()
