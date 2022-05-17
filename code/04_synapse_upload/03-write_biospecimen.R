#   Create the 'biospecimen.csv' file, 1 of 3 required metadata files for upload
#   to synapse

library('sessioninfo')
library('here')
library('readxl')

write_path = here('processed-data', '04_synapse_upload', 'biospecimen.csv')
pd_path = here('processed-data', '04_synapse_upload', 'pd.csv')
template_path = here(
    'raw-data', 'synapse_templates', 'template_biospecimen.xlsx'
)

###############################################################################
#  Populate biospecimen data frame
###############################################################################

pd = read.csv(pd_path)

meta_df = data.frame(
    'individualID' = pd$subject,
    'specimenID' = pd$Sample,
    'organ' = 'brain',
    'organWeight' = NA,
    'organRIN' = NA,
    'tissue' = 'dorsolateral prefrontal cortex',
    'isPostMortem' = TRUE,
    'BrodmannArea' = NA,
    'nucleicAcidSource' = stop(), # need to ask about this
    'cellType' = NA,
    'reprogrammedCellType' = NA,
    'terminalDifferentiationPoint' = NA,
    'passage' = NA,
    'samplingAge' = pd$age,
    'samplingAgeUnits' = 'years',
    'sampleStatus' = NA
)

#  Ensure we included all the required columns
template_names = colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names=FALSE)

session_info()
