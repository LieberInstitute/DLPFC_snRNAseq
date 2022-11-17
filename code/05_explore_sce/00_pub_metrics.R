
library("SingleCellExperiment")
library("here")

#### Load SCE ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

dim(sce)
# [1] 36601 77604

table(sce$Sample)
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 3101        5911        2861        2723        5205        4282        3898        4067        3059        3212 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 4724        4004        4275        4707        4020        4997        2661        5774        4123

summary(as.numeric(table(sce$Sample)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2661    3156    4067    4084    4716    5911 


summary(as.numeric(table(sce$Sample[sce$cellType_hc != "drop"])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 921    2279    2928    2971    3686    5043

## What nuc are dropped?
table(sce$cellType_hc != "drop")

# FALSE  TRUE 
# 21157 56447 

table(is.na(sce$cellType_layer), sce$cellType_hc != "drop")

#       FALSE  TRUE
# FALSE     0 54394
# TRUE  21157  2053

table(sce$cellType_layer)
# Astro    EndoMural        Micro        Oligo          OPC   Excit_L2/3     Excit_L3 Excit_L3/4/5     Excit_L4 
# 3979         2157         1601        10894         1940           82        10459         3043         2388 
# Excit_L5   Excit_L5/6     Excit_L6        Inhib 
# 2505         2487         1792        11067

## Total UMIs
summary(sce$sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 220    2064    4841   13660   17590  296099

## non zero-expression per nuc
zero_sums <- colSums(counts(sce) != 0)
summary(zero_sums)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 209    1403    2543    3605    5470   14808

