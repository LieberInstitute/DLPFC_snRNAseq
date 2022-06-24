## adapted from https://www.stephaniehicks.com/biocdemo/articles/Demo.html

library("SingleCellExperiment")
# library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")
# library("jaffelab")

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

# load(here("processed-data", "03_build_sce","km_res.Rdata"),verbose = TRUE)

#### Use fasthplus + wcss to refine k ####
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
fasthplus <- lapply(l, function(li) {
  message(Sys.time())
  initial_t <- find_t(L=li,proportion = 0.01)
  h <- hpb(D=reducedDims(sce)$HARMONY,
           L=li,
           t=initial_t,
           r=30)
})

fasthplus <- unlist(fasthplus)

## get wcss
wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

km_metrics <- data.frame (k=k_list, wcss = wcss, fasthplus = fasthplus)
write.csv(km_metrics, file = here("processed-data","03_build_sce","mb_kmeans_metrics.csv"))

#### Plot metrics to select best k #### 
pdf(here("plots","03_build_sce","cluster", "mb_kmeans_wcss-fastH.pdf"))
plot(k_list, wcss, type = "b")
abline(v=29, lty=2, col="red")
plot(k_list, fasthplus, type = "b")
abline(v=29, lty=2, col="red")
dev.off()

## Save data
save(km_res,km_metrics, file = here("processed-data", "03_build_sce","km_res.Rdata"))

message("Saving Data - ", Sys.time())
saveHDF5SummarizedExperiment(sce, dir = here("processed-data", "sce","sce_DLPFC"))


# sgejobs::job_single('cluster_mb_kmeans', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript cluster_mb_kmeans.R")
## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
