#  Create the 'assay.csv' file, 1 of 3 required metadata files for upload to
#  synapse

library('sessioninfo')
library('here')
library('readxl')

write_path = here('processed-data', '04_synapse_upload', 'assay.csv')
pd_path = here('processed-data', '04_synapse_upload', 'pd.csv')
template_path = here(
    'raw-data', 'synapse_templates', 'template_assay_rnaSeq.xlsx'
)

###############################################################################
#  Populate assay data frame
###############################################################################

pd = read.csv(pd_path)

meta_df = data.frame(
    'specimenID' = pd$Sample,
    'libraryID' = NA,
    'assay' = 'scrnaSeq',
    'platform' = stop(),
    'sampleBarcode' = stop(),
    'RIN' = NA,
    'referenceSet' = 'GRCh38',
    'rnaBatch' = NA,
    'libraryBatch' = NA,
    'sequencingBatch' = NA,
    'libraryPrep' = stop(),
    'libraryPreparationMethod' = NA,
    'libraryVersion' = stop(),
    'libraryType' = stop(),
    'isStranded' = stop(),
    'readStrandOrigin' = NA,
    'readLength' = stop(),
    'runType' = 'pairedEnd',
    'totalReads' = stop(),
    'validBarcodeReads' = stop(),
    'numberCells' = stop(),
    'medianGenes' = stop(),
    'medianUMIs' = stop()
)

#  Ensure we included all the required columns
template_names = colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names=FALSE)

session_info()
