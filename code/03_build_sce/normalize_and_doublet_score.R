library("SingleCellExperiment")
library("scran")
library("scater")
library("scuttle")
library("scDblFinder")
library("jaffelab")
library("tidyverse")
library("batchelor")
library("here")
library("sessioninfo")

# Plotting themes
my_theme <- theme_bw() +
  theme(text = element_text(size=15))

## Load empty-free sce data
load(here("processed-data", "sce", "sce_no_empty_droplets.Rdata"), verbose = TRUE)
dim(sce)
# [1] 36601 84756

sample_data <- colData(sce) %>%
  as_tibble() %>%
  select(Sample, subject, region, region_short, round) %>%
  unique() 

sample_data %>% count(round)

#### Check for low quality nuc ####
## maybe do in droplet_qc?

## High mito
# sce$high.mito.sample ## standard name?
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads=3, type="higher", batch = sce$Sample)
table(sce$high_mito)
# FALSE  TRUE 
# 78197  6559 
table(sce$high_mito, sce$Sample)


## low library size
# sce$qc.lib ## standard name?
sce$low_sum <- isOutlier(sce$sum, log=TRUE, type="lower", batch = sce$Sample)
table(sce$low_sum)
# FALSE  TRUE 
# 84294   462


## low detected features
# sce$qc.detected
sce$low_detected <- isOutlier(sce$detected, log=TRUE, type="lower", batch = sce$Sample)
table(sce$low_detected)
# FALSE  TRUE 
# 83513  1243

## All low sum are also low detected
table(sce$low_sum, sce$low_detected)
#       FALSE  TRUE
# FALSE 83513   781
# TRUE      0   462

## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_sum | sce$low_detected

table(sce$discard_auto)
# FALSE  TRUE 
# 77604  7152

## discard 8% of nuc
100*sum(sce$discard_auto)/ncol(sce)
# [1] 8.438341

(qc_t <- addmargins(table(sce$Sample, sce$discard_auto)))
#             FALSE  TRUE   Sum
# Br2720_mid   3101   642  3743
# Br2720_post  5911    39  5950
# Br2743_ant   2861   337  3198
# Br2743_mid   2723   703  3426
# Br3942_ant   5205  1064  6269
# Br3942_mid   4282    27  4309
# Br6423_ant   3898   101  3999
# Br6423_post  4067     0  4067
# Br6432_ant   3059   373  3432
# Br6471_ant   3212   410  3622
# Br6471_mid   4724   628  5352
# Br6522_mid   4004    31  4035
# Br6522_post  4275    20  4295
# Br8325_ant   4707  1244  5951
# Br8325_mid   4020   970  4990
# Br8492_mid   4997   192  5189
# Br8492_post  2661   328  2989
# Br8667_ant   5774    32  5806
# Br8667_mid   4123    11  4134
# Sum         77604  7152 84756

round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#             FALSE  TRUE   Sum
# Br2720_mid   82.8  17.2 100.0
# Br2720_post  99.3   0.7 100.0
# Br2743_ant   89.5  10.5 100.0
# Br2743_mid   79.5  20.5 100.0
# Br3942_ant   83.0  17.0 100.0
# Br3942_mid   99.4   0.6 100.0
# Br6423_ant   97.5   2.5 100.0
# Br6423_post 100.0   0.0 100.0
# Br6432_ant   89.1  10.9 100.0
# Br6471_ant   88.7  11.3 100.0
# Br6471_mid   88.3  11.7 100.0
# Br6522_mid   99.2   0.8 100.0
# Br6522_post  99.5   0.5 100.0
# Br8325_ant   79.1  20.9 100.0
# Br8325_mid   80.6  19.4 100.0
# Br8492_mid   96.3   3.7 100.0
# Br8492_post  89.0  11.0 100.0
# Br8667_ant   99.4   0.6 100.0
# Br8667_mid   99.7   0.3 100.0
# Sum          91.6   8.4 100.0

#### QC plots ####
pdf(here("plots","03_build_sce","QC_outliers.pdf"), width = 21)
## Mito rate
plotColData(sce, x = "Sample", y="subsets_Mito_percent", colour_by="high_mito") +
  ggtitle("Mito Precent")+
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)

# ## low sum
plotColData(sce, x = "Sample", y="sum", colour_by="low_sum") +
  scale_y_log10() +
  ggtitle("Total count") +
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

# ## low detected
plotColData(sce, x = "Sample", y="detected", colour_by="low_detected") +
  scale_y_log10() +
  ggtitle("Detected features") +
  # geom_hline(yintercept = 500)+
  facet_wrap(~sce$round, scales = "free_x", nrow = 1)

# Mito rate vs n detected features

  plotColData(sce, x="detected", y="subsets_Mito_percent",
              colour_by="discard_auto", point_size=2.5, point_alpha=0.5) 

# Detected features vs total count

  plotColData(sce, x="sum", y="detected",
              colour_by="discard_auto", point_size=2.5, point_alpha=0.5) 

dev.off()


#### Doublet detection / removal ####
## To speed up, run on sample-level top-HVGs - just take top 1000 
set.seed(506)

dbl_dens_sample <- map(splitit(sce$Sample), function(s){
  sce_temp <- sce[,s]
  ## To speed up, run on sample-level top-HVGs - just take top 1000 
  normd <- logNormCounts(sce_temp)
  geneVar <- modelGeneVar(normd)
  topHVGs <- getTopHVGs(geneVar, n=1000)
  
  dbl_dens <- computeDoubletDensity(normd, subset.row=topHVGs)
  return(dbl_dens)
})

map_dbl(dbl_dens_sample, median)
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 0.314412    1.047200    0.377364    0.267228    0.376140    0.448136    0.351912    0.300958    0.480480    0.347712 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 0.278304    1.654350    1.812490    0.297550    0.249500    0.456632    0.514108    0.429644    0.471276 

map(dbl_dens_sample, ~round(quantile(.x, probs=seq(0,1,by=0.05)),1))

## Visualize doublet scores ##

dbl_df <- stack(dbl_dens_sample) %>%
  rename(dbl_score = values, Sample = ind)

dbl_box_plot <- dbl_df %>%
  ggplot(aes(x = reorder(Sample, dbl_score, FUN = median), y = dbl_score)) +
  geom_boxplot() +
  labs(x = "Sample") +
  geom_hline(yintercept = 5, color = "red", linetype = "dashed") +
  coord_flip() +
  my_theme

ggsave(dbl_box_plot, filename = here("plots","03_build_sce","doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
  ggplot(aes(x  = dbl_score)) +
  geom_density() +
  labs(x = "doublet score") +
  facet_grid(Sample ~ .) +
  theme_bw()

ggsave(dbl_density_plot, filename = here("plots","03_build_sce","doublet_scores_desnity.png"), height = 17)


# Percent that would be dropped at density score >= 5
map_int(dbl_dens_sample,  ~sum(.x >= 5))
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 45          97          21          11          52           5          30          23          32          37 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 49         135         366          60          31          10           4          39          12 
map_dbl(dbl_dens_sample,  ~round(100*sum(.x >= 5)/length(.x), 1))
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 1.2         1.6         0.7         0.3         0.8         0.1         0.8         0.6         0.9         1.0 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 0.9         3.3         8.5         1.0         0.6         0.2         0.1         0.7         0.3 


cn_test <- map(splitit(sce$Sample), function(s){
  sce_temp <- sce[,s]
  return(colnames(sce_temp))})

cn_test2 <- unlist(cn_test)[unlist(splitit(sce$Sample))]

table(cn_test2 == colnames(sce))
table(head(unlist(cn_test)) == head(colnames(sce)))
all(unlist(cn_test) == colnames(sce))

