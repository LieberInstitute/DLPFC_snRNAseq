BiocManager::install("lahuuki/DeconvoBuddies")

library("SingleCellExperiment")
library("scater")
library("here")
library("sessioninfo")
library("DeconvoBuddies")

load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)

#### Find Markers ####
table(sce$cellType_hc)
# Astro Endo.Mural_01 Endo.Mural_02      Excit_01      Excit_02      Excit_03      Excit_04      Excit_05 
# 3979           446          1711          7927          2487          1309          2171          2532 
# Excit_06      Excit_07      Excit_08      Excit_09      Excit_10      Excit_11      Excit_12      Excit_13 
# 329           334          1463          2561          1079           482           420          1567 
# Excit_14      Excit_15      Inhib_01      Inhib_02      Inhib_03      Inhib_04      Inhib_05      Inhib_06 
# 82            66          5366          1267          1310           565          1192          1367 
# Micro      Oligo_01      Oligo_02      Oligo_03           OPC 
# 1601         23025          4732          4294          1940 

colLabels(sce) <- sce$cellType_hc

table(sce$Sample)
mod <- "~Sample"

# message("Running 1vALL findMarkers")
# Sys.time()
# 
# markers_1vALL_sample <- findMarkers_1vAll(sce, assay_name = "logcounts", cellType_col = "cellType_hc", mod = mod)
# 
# message("Done - ", Sys.time())
# Sys.time()
# 
# save(markers_1vALL_sample, file = here("processed-data", "03_build_sce","cell_type_markers_1vALL_mod.Rdata"))

## subtype enrichment 

## Which cell types have sub-types?
cell_types <- levels(sce$cellType_hc)
(ctt <- table(jaffelab::ss(cell_types,"_")))
# Astro EndoMural     Excit     Inhib     Micro     Oligo       OPC 
# 1         2        15         6         1         3         1

## 
(sub_cell_types <- names(ctt)[ctt > 1])
# [1] "EndoMural" "Excit"     "Inhib"     "Oligo" 

markers_enrich_subtype <- purrr::map(sub_cell_types, function(ct){
  
  message("Finding ", ct, " enrichment markers - ", Sys.time())
  
  sce_temp <- sce[,grepl(ct, sce$cellType_hc)]
  sce_temp$cellType_hc <- droplevels(sce_temp$cellType_hc)
  
  print(table(sce_temp$cellType_hc))
  enrich_subtype <- findMarkers_1vAll(sce_temp, assay_name = "logcounts", cellType_col = "cellType_hc", mod = mod)
  
  message("Done - ", Sys.time())
  
  return(enrich_subtype)
})

save(markers_enrich_subtyp, file = here("processed-data", "03_build_sce","cell_type_markers_enrich_subtype.Rdata"))


# sgejobs::job_single('09_find_markers2', create_shell = TRUE, queue= 'bluejay', memory = '20G', command = "Rscript 09_find_markers2.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

