library(tidyverse)
library(SingleCellExperiment)


# Receive Arguments from Command Line
args=(commandArgs(TRUE))

if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        print(args[[i]])
        eval(parse(text=args[[i]]))
    }
}

## NOTE: reserved for when command line arguments fails
# crn_Positions <- c("Anterior", "Middle", "Posterior")
# crn_sec <- crn_Positions[1]


# Create section specific folders to contain results
# TODO: change this to LIBD drive
if(!dir.exists("~/CCC_snRNA")) dir.create("~/CCC_snRNA")
fdl_path <- paste0("~/CCC_snRNA/", crn_sec, "/")
dir.create(fdl_path, recursive = TRUE)

## NOTE: test
print(crn_sec)



# Data Prep ---------------------------------------------------------------

# Load Proprocessed Data
## NOTE: Set working directory is necessary in order to use the rds data because of the .h5 file
## https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1668537164053579?thread_ts=1668530634.509189&cid=C01EA7VDJNT
# setwd("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/sce/sce_DLPFC_annotated/")
# sce <- readRDS("se.rds")
#
# # Sample Summary Statistics
# # col_data <- colData(sce) |> as.data.frame()
# # col_data |> group_by(BrNum) |> summarize(n = n_distinct(pos))
#
# # Subset cells from each coronal section
# sce_crn <- sce[, sce$Position==crn_sec]
#
# # if(ncol(sce_crn) <=0) stop
#
# # QC: Remove bad cells----------------------------------------------------------------------
# # TODO: ask Louise why
# # https://jhu-genomics.slack.com/archives/C01EA7VDJNT/p1668530673111519?thread_ts=1668530634.509189&cid=C01EA7VDJNT
# sce_crn <- sce_crn[,sce_crn$cellType_hc != "drop"]
# sce_crn$cellType_hc <- droplevels(sce_crn$cellType_hc)
#
# ## See if batch effect is large
# # library(scater)
# # ggsave(paste0( fdl_path,"tsne.png"),
# #        plot = plotReducedDim(sce_crn, dimred="TSNE", colour_by="BrNum"))
#
# ## Empty Cells/Genes
# # sum(rowSums(logcounts(sce_crn))==0)
# # sum(colSums(logcounts(sce_crn))==0)
#
#
# # HVG
# # TODO: Decides if we want to prioritize some genes to reduce computation load
# # TODO: If so, change the sce object names
# # library(scran)
# # dec_sce <- modelGeneVar(sce_Br2720)
# # hvg_names <- getTopHVGs(dec_sce, n=100)
#
# # LIANA Analysis ----------------------------------------------------------
# library(liana)
# # Prep
# ## Assigning cellType to sce label (Required by LIANA)
# colLabels(sce_crn) <- colData(sce_crn)$cellType_broad_hc
#
# # TEST: small scale
# # TODO: remove this when scaling up the analysis
# sce_crn <- sce_crn[, 1:1000]
#
#
# # Run liana
# liana_test <- liana_wrap(sce_crn)
# # TODO: save this result
#
#
# # If segfault happens, remove that list in the list.
# #  tmp <- liana_test
# # tmp$cellphonedb <- NULL
#
# # Consensus Ranks
# # TODO: save this result
# liana_res <- liana_test %>%
#     liana_aggregate()
#
#
# # Plot
# # library(circlize)
# # pdf(file = "~/CCC_toy.pdf")
# # chord_freq(tmp#,
# #                 # source_groups = c("CD8 T", "NK", "B"),
# #                 # target_groups = c("CD8 T", "NK", "B")
# #                 )
# # dev.off()
