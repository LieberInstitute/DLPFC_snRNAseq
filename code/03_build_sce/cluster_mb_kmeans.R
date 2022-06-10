## adapted from https://www.stephaniehicks.com/biocdemo/articles/Demo.html

library("SingleCellExperiment")
library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")
library("numform")

## load data
load(here("processed-data", "03_build_sce","sce_harmony_Sample.Rdata"), verbose = TRUE)

## neron paper returned 19 clusters for DLPFC, try 5:50
set.seed(610)
k_list <- seq(5, 50)

km_res <- lapply(k_list, function(k) {
  message("k=",k)
  mbkmeans(sce, clusters = k, 
           batch_size = 500,
           reduceMethod = "HARMONY",
           calc_wcss = TRUE)
})

names(km_res[[1]])

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

fasthplus <- lapply(l, function(li) {
  message(Sys.time())
  initial_t <- find_t(L=li,proportion = 0.01)
  h <- hpb(D=reducedDims(sce)$HARMONY,
           L=li,
           t=initial_t,
           r=30)
})

wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))
results <- data.frame (k=k_list, wcss = wcss, fasthplus = fasthplus)
write.table(results,file = here::here("processed-data","rdata","spe","06_fasthplus","fasthplus_results_no_WM.csv"), append = TRUE)

pdf(here("plots","03_build_sce","cluster", "mb_kmeans_wcss.pdf"))
plot(k_list, wcss, type = "b")
abline(v=26, lty=2, col="red")
dev.off()

sce$kmeans_26 <- paste0("mbk", f_pad_zero(km_res[[which(k_list==26)]]$Clusters))
table(sce$kmeans_26)
# mbk01 mbk02 mbk03 mbk04 mbk05 mbk06 mbk07 mbk08 mbk09 mbk10 mbk11 mbk12 mbk13 mbk14 mbk15 mbk16 mbk17 mbk18 mbk19 mbk20 
# 1690  2127  3777  3800  1686     3  1280 20016   533    14  2144  4399  1516  2859   304  2759  1081  3829   608 13377 
# mbk21 mbk22 mbk23 mbk24 mbk25 mbk26 
# 1294  4729     2   432  1968  1377 

sce$kmeans_30 <- paste0("mbk", f_pad_zero(km_res[[which(k_list==30)]]$Clusters))
table(sce$kmeans_30)

# mbk01 mbk02 mbk03 mbk04 mbk05 mbk06 mbk07 mbk08 mbk09 mbk10 mbk11 mbk12 mbk13 mbk14 mbk15 mbk16 mbk17 mbk18 mbk19 mbk20 
# 366  5740  2571  4707    26  2895   963   381  1738  2466  1589    17  2895   289   217     2   404   881 23842 12300 
# mbk21 mbk22 mbk23 mbk24 mbk25 mbk26 mbk27 mbk28 mbk29 mbk30 
# 1243  1408  5030   705    29  1018  1239  1098   490  1055 


## Plotting 

# Plotting set up 
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

plot_dir = here("plots","03_build_sce","cluster")


clusters_k26 <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=kmeans_26)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  my_theme + 
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(clusters_k26, filename = here(plot_dir, "clusters_mbkm-26.png"), width = 10)


clusters_k30 <- ggcells(sce, mapping=aes(x=TSNE.1, y=TSNE.2, colour=kmeans_30)) +
  geom_point(size = 0.2, alpha = 0.3) +
  coord_equal() +
  my_theme + 
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave(clusters_k30, filename = here(plot_dir, "clusters_mbkm-30.png"), width = 10)
  
