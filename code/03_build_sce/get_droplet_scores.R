## Based on
## https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/10x_all-FACS-n10_2021rev_step01_processing-QC_MNT.R

library("SingleCellExperiment")
library("DropletUtils")
library("BiocParallel")
library("scuttle")
library("jaffelab")
library("here")
library("sessioninfo")

## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)

table(sce$Sample)
# Br2720_mid Br2720_post  Br2743_ant  Br2743_mid  Br3942_ant  Br3942_mid  Br6423_ant Br6423_post  Br6432_ant  Br6471_ant 
# 1152093      922885      655424     1495673     1692299     2209670     1521781     1786747      736110     1148322 
# Br6471_mid  Br6522_mid Br6522_post  Br8325_ant  Br8325_mid  Br8492_mid Br8492_post  Br8667_ant  Br8667_mid 
# 1653348      604278      681586     1855444     2192531     1302254      965917     1875106     2154338

sample_index <- splitit(sce$Sample)

message("Get Droplett scores by sample")

Sys.time()
e.out <- lapply(sample_index[c(1,2)], function(i){
  sce_sample <- sce[, i]
  e.sample <- DropletUtils::emptyDrops(
    sce,
    niters = 30000
    # ,
    # BPPARAM = BiocParallel::MulticoreParam(4)
  )
  return(e.sample)
})
Sys.time()

save(e.out, file = here("processed-data", "03_build_sce", "droplet_scores.Rdata"))

# sgejobs::job_single('get_droplet_scores', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript get_droplet_scores.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
