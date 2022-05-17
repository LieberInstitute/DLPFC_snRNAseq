#  Create the 'manifest.tsv' file, the "main" file required for upload to
#  synapse

library('sessioninfo')
library('here')
library('readxl')

write_path = here('processed-data', '04_synapse_upload', 'manifest.tsv')
pd_path = here('processed-data', '04_synapse_upload', 'pd.csv')
template_path = here(
    'raw-data', 'synapse_templates', 'template_manifest.xlsx'
)
fastq_naming_path = here(
    'processed-data', '04_synapse_upload', 'fastq_renaming_scheme.csv'
)

pd = read.csv(pd_path)

paths = c( 
    here('processed-data', '04_synapse_upload', 'assay.csv'),
    here('processed-data', '04_synapse_upload', 'biospecimen.csv'),
    here('processed-data', '04_synapse_upload', 'individual.csv'),
    fastq_naming$new_path
)

###############################################################################
#  Populate data frame that will become TSV
###############################################################################

num_metadata = 3
num_fastq = length(fastq_naming$new_path)

#  Populate a data frame
meta_df = data.frame(
    'path' = paths,
    'parent' = stop(), # need to ask this
    'individualID' = pd$subject[
        match(fastq_naming$sample_id, pd$Sample)
    ],
    'specimenID' = fastq_naming$sample_id,
    'isMultiIndividual' = c(
        rep(NA, num_metadata), rep(FALSE, num_fastq)
    ),
    'isMultiSpecimen' = c(
        rep(NA, num_metadata), rep(FALSE, num_fastq)
    ),
    'assay' = 'scrnaSeq',
    'libraryID' = libraryID,
    'fileFormat' = c(
        rep('csv', num_metadata), rep('fastq', num_fastq)
    ),
    'consortium' = 'PEC',
    'study' = stop(), # need to ask this
    'grant' = stop(), # need to ask this
    'resourceType' = 'experimentalData',
    'dataType' = NA,
    'dataSubtype' = c(
        rep('metadata', num_metadata), rep('raw', num_fastq)
    ),
    'metadataType' = metadataType = c(
        'assay', 'biospecimen', 'individual',
        rep(NA, num_fastq)
    ),
    'analysisType' = NA
)

#   Ensure we included all the required columns
template_names = colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.table(meta_df, write_path, row.names=FALSE, sep = '\t')

session_info()
