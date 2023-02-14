
library("SingleCellExperiment")
library("here")
library("tidyverse")

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

#### Supp Table export ####

list.files(here("processed-data","03_build_sce"))
list.files(here("processed-data","02_cellranger_metrics"))

sample_info <- read.csv(here("processed-data","03_build_sce","sample_info.csv")) |>
  rename(Donor = subject, Position = region)

droplet_summary <- read.csv(here("processed-data","03_build_sce","drop_summary.csv")) |>
  select(Sample, total_droplets = total_n, droplet_lower_cutoff = lower_cutoff, pre_QC_n = non_empty)

cell_ranger_info <-read.csv(here("processed-data","02_cellranger_metrics","cellranger_metrics.csv")) |>
  select(file_id = X, Number.of.Reads, Mean.Reads.per.Cell, Median.UMI.Counts.per.Cell, Median.Genes.per.Cell) #mean unique molecular indices - only median in this table

colnames(sample_info_all)
# [1] "Sample"                                         "file_id"                                       
# [3] "region"                                         "subject"                                       
# [5] "round"                                          "n"                                             
# [7] "n_high_mito"                                    "n_low_sum"                                     
# [9] "n_low_detect"                                   "n_discard_auto"                                
# [11] "Estimated.Number.of.Cells"                      "Mean.Reads.per.Cell"                           
# [13] "Median.Genes.per.Cell"                          "Number.of.Reads"                               
# [15] "Valid.Barcodes"                                 "Sequencing.Saturation"                         
# [17] "Q30.Bases.in.Barcode"                           "Q30.Bases.in.RNA.Read"                         
# [19] "Q30.Bases.in.UMI"                               "Reads.Mapped.to.Genome"                        
# [21] "Reads.Mapped.Confidently.to.Genome"             "Reads.Mapped.Confidently.to.Intergenic.Regions"
# [23] "Reads.Mapped.Confidently.to.Intronic.Regions"   "Reads.Mapped.Confidently.to.Exonic.Regions"    
# [25] "Reads.Mapped.Confidently.to.Transcriptome"      "Reads.Mapped.Antisense.to.Gene"                
# [27] "Fraction.Reads.in.Cells"                        "Total.Genes.Detected"                          
# [29] "Median.UMI.Counts.per.Cell"

## from manuscript
# droplets were sequenced


sample_info_all <- sample_info |> 
  left_join(droplet_summary) |>
  left_join(cell_ranger_info) |>
  select(Sample, file_id, Donor, Position, round, ## sample info
         Number.of.Reads, Mean.Reads.per.Cell, Median.UMI.Counts.per.Cell, Median.Genes.per.Cell, ## cell ranger
         total_droplets, droplet_lower_cutoff, pre_QC_n,  ## empty drops 
         n_high_mito, n_low_sum, n_low_detect, n_discard_auto) |> ## QC
  mutate(n = pre_QC_n - n_discard_auto)

sample_info_all |> 
  select(Median.Genes.per.Cell)|>
  summary()

# Sample            file_id             Donor             Position            round           Number.of.Reads    
# Length:19          Length:19          Length:19          Length:19          Length:19          Min.   :162821623  
# Class :character   Class :character   Class :character   Class :character   Class :character   1st Qu.:275981256  
# Mode  :character   Mode  :character   Mode  :character   Mode  :character   Mode  :character   Median :323114123  
#                                                                                                Mean   :351032590  
#                                                                                                3rd Qu.:483252608  
#                                                                                                Max.   :541289607  
# Mean.Reads.per.Cell Median.UMI.Counts.per.Cell Median.Genes.per.Cell total_droplets    droplet_lower_cutoff
# Min.   : 38062      Min.   : 1323              Min.   : 985          Min.   : 604278   Min.   :219.0       
# 1st Qu.: 47784      1st Qu.: 2144              1st Qu.:1434          1st Qu.: 944401   1st Qu.:276.5       
# Median : 55811      Median : 3334              Median :1959          Median :1495673   Median :390.0       
# Mean   : 65253      Mean   : 5253              Mean   :2415          Mean   :1400306   Mean   :441.3       
# 3rd Qu.: 75312      3rd Qu.: 7366              3rd Qu.:3259          3rd Qu.:1821096   3rd Qu.:581.0       
# Max.   :136799      Max.   :14627              Max.   :5146          Max.   :2209670   Max.   :910.0       
# pre_QC_n     n_high_mito       n_low_sum       n_low_detect    n_discard_auto         n       
# Min.   :2989   Min.   :   0.0   Min.   :  0.00   Min.   :  0.00   Min.   :   0.0   Min.   :2661  
# 1st Qu.:3682   1st Qu.:  26.5   1st Qu.:  0.00   1st Qu.:  0.00   1st Qu.:  31.5   1st Qu.:3156  
# Median :4134   Median : 258.0   Median :  0.00   Median :  0.00   Median : 328.0   Median :4067  
# Mean   :4461   Mean   : 345.2   Mean   : 24.32   Mean   : 65.42   Mean   : 376.4   Mean   :4084  
# 3rd Qu.:5270   3rd Qu.: 565.0   3rd Qu.:  0.00   3rd Qu.: 70.50   3rd Qu.: 635.0   3rd Qu.:4716  
# Max.   :6269   Max.   :1229.0   Max.   :192.00   Max.   :412.00   Max.   :1244.0   Max.   :5911

library(xlsx)

key <- data.frame(
  colnames = colnames(sample_info_all)
)

## Clear file and write key
annotation_xlsx <-
  here("processed-data","05_explore_sce", "00_pub_metrics", "sn_sample_info.xlsx")

write.xlsx(
  key,
  file = annotation_xlsx,
  sheetName = "Key",
  append = FALSE,
  row.names = FALSE
)

## write annotations
write.xlsx(
  sample_info_all,
  file = annotation_xlsx,
  sheetName = paste0("snRNA-seq_sample_info"),
  append = TRUE,
  row.names = FALSE
)
