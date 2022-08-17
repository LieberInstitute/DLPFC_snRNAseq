library("SingleCellExperiment")
library("scran")
library("here")

#### Load data ####
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

pd <- as.data.frame(SummarizedExperiment::colData(sce))
# message("donor" %in% colnames(pd))

mod <- with(pd, stats::model.matrix(~Sample))
mod <- mod[, -1, drop = F] # intercept otherwise automatically dropped by `findMarkers()`
head(mod)

sce$contrast <- ifelse(sce[["cellType_hc"]] == "Astro", 1, 0)
table(sce$contrast)
# 0     1 
# 73625  3979 

fm_astro <- scran::findMarkers(sce,
                         groups = sce$contrast,
                         assay.type = "logcounts", design = mod, test.type = "t",
                         direction = "up", pval.type = "all", full.stats = T
)

save(fm_astro, file =  here("processed-data","03_build_sce","10_markere_deep_dive","fm_astro.Rdata"))

##  sorted marker gene list for the corresponding group -> so cellType.target == 1 ...? 
fm_astro[["0"]]
# DataFrame with 6 rows and 4 columns
# p.value       FDR summary.stats                   stats.1
# <numeric> <numeric>     <numeric>               <DataFrame>
# NRXN3            0         0       1.77907 1.77907:-3915.14:-3904.63
# ANK3             0         0       1.16963 1.16963:-3706.92:-3697.11
# PTPRD            0         0       1.55819 1.55819:-3678.14:-3668.73
# CNTNAP2          0         0       1.56137 1.56137:-3420.35:-3411.23
# IL1RAPL1         0         0       1.47912 1.47912:-3259.62:-3250.72
# RBFOX1           0         0       1.67179 1.67179:-3038.95:-3030.24

## Returned by my code
head(fm_astro[["1"]])
# DataFrame with 6 rows and 4 columns
# p.value       FDR summary.stats                   stats.0
# <numeric> <numeric>     <numeric>               <DataFrame>
# ADGRV1             0         0       2.13884 2.13884:-29283.4:-29272.9
# OBI1-AS1           0         0       1.67705 1.67705:-26898.8:-26889.0
# LINC00299          0         0       1.27898 1.27898:-20584.9:-20575.5
# AL137139.2         0         0       1.09167 1.09167:-20459.0:-20449.9
# SLC1A2             0         0       2.75541 2.75541:-19783.7:-19774.8
# PITPNC1            0         0       2.04256 2.04256:-19769.0:-19760.3

summary(fm_astro[[1]])

fm_astro[["1"]]$stats.0
# DataFrame with 36601 rows and 3 columns
# logFC log.p.value   log.FDR
# <numeric>   <numeric> <numeric>
# ADGRV1      2.13884    -29283.4  -29272.9 ## found in a paper as ASTRO marker
# OBI1-AS1     1.67705    -26898.8  -26889.0
# LINC00299    1.27898    -20584.9  -20575.5
# AL137139.2   1.09167    -20459.0  -20449.9
# SLC1A2       2.75541    -19783.7  -19774.8 ## marker we use for annotation
# ...              ...         ...       ...
# PAK3       -0.432842           0         0
# GRIA3      -0.487904           0         0
# AL008633.1 -0.402078           0         0
# LINC00632  -0.478628           0         0
# IDS        -0.454775           0         0

# std.lfc=TRUE, the log-fold change for each gene is standardized by the variance
# Standardized log-fold changes may be more appealing for visualization as it avoids large fold changes due to large variance. The choice
# of std.lfc does not affect the calculation of the p-values. aka t-stat
fm_std_astro <- scran::findMarkers(sce,
                             groups = sce$contrast,
                             assay.type = assay_name, design = mod, test.type = "t",
                             std.lfc = TRUE,
                             direction = "up", pval.type = "all", full.stats = T
)


pw_astro <- findMarkers(sce, groups=sce$cellType_hc,
                            assay.type="logcounts", design=mod, test="t",
                            direction="up", pval.type="all", full.stats=T)

save(pw_astro, file =  here("processed-data","03_build_sce","10_markere_deep_dive","pw_astro.Rdata"))
