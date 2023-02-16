library("SingleCellExperiment")
library("scater")
library("here")
library("sessioninfo")
library("DeconvoBuddies")
library("purrr")

load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## Which cell types have sub-types?
cell_types <- levels(sce$cellType_hc)
(ctt <- table(jaffelab::ss(cell_types, "_")))
# Astro EndoMural     Excit     Inhib     Micro     Oligo       OPC
# 1         2        15         6         1         3         1

##
(sub_cell_types <- names(ctt)[ctt > 1])
# [1] "EndoMural" "Excit"     "Inhib"     "Oligo"

names(sub_cell_types) <- sub_cell_types

colLabels(sce) <- sce$cellType_hc

subtype_markers <- map(sub_cell_types, function(ct) {
    message("Finding ", ct, " markers - ", Sys.time())

    sce_temp <- sce[, grepl(ct, sce$cellType_hc)]
    sce_temp$cellType_hc <- droplevels(sce_temp$cellType_hc)

    print(table(sce_temp$cellType_hc))
    markers_all <- scran::findMarkers(sce_temp, pval.type = "all", direction = "up")
    markers_any <- scran::findMarkers(sce_temp, pval.type = "any", direction = "up")

    message("Done - ", Sys.time())

    return(list(all = markers_all, any = markers_any))
})


save(subtype_markers, file = here("processed-data", "03_build_sce", "subtype_markers.Rdata"))

# sgejobs::job_single('11_advanced_marker_search', create_shell = TRUE, memory = '100G', command = "Rscript 11_advanced_marker_search.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
