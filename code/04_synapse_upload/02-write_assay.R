#  Create the 'assay.csv' file, 1 of 3 required metadata files for upload to
#  synapse. We need this type of file for the scRNA-seq FASTQs and another for
#  the two genotyping files.

library("sessioninfo")
library("here")
library("readxl")
library("jaffelab")

write_rna_path <- here(
    "processed-data", "04_synapse_upload", "assay_scrnaSeq.csv"
)
write_snp_path <- here(
    "processed-data", "04_synapse_upload", "assay_snpArray.csv"
)
pd_path <- here("processed-data", "04_synapse_upload", "pd.csv")
template_rna_path <- here(
    "raw-data", "synapse_templates", "template_assay_rnaSeq.xlsx"
)
template_snp_path <- here(
    "raw-data", "synapse_templates", "template_assay_snpArray.xlsx"
)
geno_path <- here(
    "raw-data", "sample_info", "Genotyping.Batch.info.22.09.21.corrected.csv"
)

################################################################################
#  Populate scRNA-seq assay data frame
################################################################################

pd <- read.csv(pd_path)

meta_df <- data.frame(
    "specimenID" = pd$Sample,
    "libraryID" = NA,
    "assay" = "scrnaSeq",
    "platform" = "IlluminaNovaseq6000",
    "sampleBarcode" = NA,
    "RIN" = NA,
    "referenceSet" = "GRCh38",
    "rnaBatch" = NA,
    "libraryBatch" = NA,
    "sequencingBatch" = NA,
    "libraryPrep" = "polyAselection",
    "libraryPreparationMethod" = "10x",
    "libraryVersion" = NA,
    "libraryType" = NA,
    "isStranded" = NA,
    "isResequenced" = NA,
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
template_names <- colnames(read_excel(template_rna_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_rna_path, row.names = FALSE)

################################################################################
#  Populate snpArray (genotyping) assay data frame
################################################################################

#   Read in genotype batch info, take just the date, and line up with pd
geno <- read.csv(geno_path)
geno$GenotypingBatch <- ss(geno$GenotypingBatch, "-", 3)
geno <- geno[match(pd$subject, geno$BrNum), ]

meta_df <- data.frame(
    "specimenID" = pd$Sample,
    "assay" = "snpArray",
    "platform" = "Illumina_Omni2pt5M",
    "dnaBatch" = NA,
    "arrayBatch" = geno$GenotypingBatch,
    "ratio260over280" = NA,
    "ratio260over230" = NA,
    "GQN" = NA
)

#  Ensure we included all the required columns
template_names <- colnames(read_excel(template_snp_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_snp_path, row.names = FALSE)

session_info()
