#   Create the 'individual.csv' file, 1 of 3 required metadata files for upload
#   to synapse

library("sessioninfo")
library("here")
library("readxl")

write_path <- here("processed-data", "04_synapse_upload", "individual.csv")
pd_path <- here("processed-data", "04_synapse_upload", "pd.csv")
template_path <- here(
    "raw-data", "synapse_templates", "template_individual_human.xlsx"
)

###############################################################################
#  Populate "individual" data frame
###############################################################################

#   Read in and tweak phenotype data to fit Synapse conventions
pd <- read.csv(pd_path)
pd$sex[pd$sex == "M"] <- "male"
pd$sex[pd$sex == "F"] <- "female"
pd$diagnosis[pd$diagnosis == "Control"] <- "control"

meta_df <- data.frame(
    "individualID" = pd$subject,
    "individualIdSource" = "LIBD",
    "species" = "Human",
    "reportedGender" = pd$sex,
    "sexChromosome" = NA,
    "race" = NA,
    "ethnicity" = NA,
    "genotypeInferredAncestry" = NA,
    "familialRelationship" = NA,
    "IQ" = NA,
    "BMI" = NA,
    "primaryDiagnosis" = pd$diagnosis,
    "primaryDiagnosisDetail" = NA,
    "otherDiagnosis" = NA,
    "otherDiagnosisDetail" = NA,
    "familyHistory" = NA,
    "familyHistoryDetails" = NA,
    "ageOnset" = NA,
    "neuropathDescription" = NA,
    "dementia" = NA,
    "CDR" = NA,
    "Braak" = NA,
    "otherMedicalDx" = NA,
    "otherMedicalDetail" = NA,
    "ageDeath" = pd$age,
    "ageDeathUnits" = "years",
    "causeDeath" = NA,
    "mannerDeath" = NA,
    "postmortemTox" = NA,
    "postmortemToxDetails" = NA,
    "postmortemToxSource" = NA,
    "medRecordTox" = NA,
    "PMICertain" = NA,
    "PMI" = NA,
    "pH" = NA
)

#  Ensure we included all the required columns
template_names <- colnames(read_excel(template_path))
stopifnot(all(template_names %in% colnames(meta_df)))
stopifnot(all(colnames(meta_df) %in% template_names))

write.csv(meta_df, write_path, row.names = FALSE)

session_info()
