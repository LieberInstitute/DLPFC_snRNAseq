#  Create the 'assay.csv' file, 1 of 3 required metadata files for upload to
#  synapse

library("sessioninfo")
library("here")
library("readxl")

write_path <- here("processed-data", "04_synapse_upload", "assay.csv")
pd_path <- here("processed-data", "04_synapse_upload", "pd.csv")
template_path <- here(
    "raw-data", "synapse_templates", "template_assay_rnaSeq.xlsx"
)

###############################################################################
#  Populate assay data frame
###############################################################################

pd <- read.csv(pd_path)

meta_df <- data.frame(
    "specimenID" = pd$Sample,
    "libraryID" = NA,
    "assay" = "scrnaSeq",
    "platform" = NA,
    "sampleBarcode" = NA,
    "RIN" = NA,
    "referenceSet" = "GRCh38",
    "rnaBatch" = NA,
    "libraryBatch" = NA,
    "sequencingBatch" = NA,
    "libraryPrep" = NA,
    "libraryPreparationMethod" = NA,
    "libraryVersion" = NA,
    "libraryType" = NA,
    "isStranded" = NA,
    "readStrandOrigin" = NA,
    "readLength" = NA,
    "runType" = "pairedEnd",
    "totalReads" = NA,
    "validBarcodeReads" = NA,
    "numberCells" = NA,
    "medianGenes" = NA,
    "medianUMIs" = NA
)

#  Ensure we included all the required columns
template_names <- colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names = FALSE)

session_info()
