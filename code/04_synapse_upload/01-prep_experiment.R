#   1. Create symbolic links to FASTQs with a consistent naming convention
#   2. Create a standalone CSV containing sample IDs and FASTQ file paths
#      for all 19 samples
#   3. Create a standalone CSV containing phenotype data with one row per
#      sample

library('sessioninfo')
library('here')
library('SingleCellExperiment')
library('jaffelab')

sce_path = here("processed-data", "sce", "sce_raw.Rdata")
pd_path = here('processed-data', '04_synapse_upload', 'pd.csv')
df_out_path = here(
    'processed-data', '04_synapse_upload', 'fastq_renaming_scheme.csv'
)
dest_dir = here('raw-data', 'FASTQ_renamed')

dir.create(dirname(pd_path), showWarnings = FALSE)
dir.create(dest_dir, showWarnings = FALSE)

load(sce_path, verbose = TRUE)

#   Take one row per sample of the phenotype data
pd = colData(sce)
pd = pd[
    match(unique(pd$Sample), pd$Sample),
    - match(c('Barcode', 'key'), colnames(pd))
]
write.csv(pd, file = pd_path, quote=FALSE, row.names = FALSE)

fastq_old = lapply(
    here('raw-data', 'FASTQ', pd$file_id), list.files, full.names = TRUE
)

#   The new naming convention for the FASTQ files to upload to synapse will be:
#
#       [sample ID]_[MATE]_[INDEX].fastq.gz
#
#   where [INDEX] is assigned randomly. For example, if there are 3 pairs of
#   FASTQ files for a given sample, one of the pairs will be randomly given the
#   index 1, and so on for index 2 and 3.

fastq_new = lapply(
    1:nrow(pd),
    function(i) {
        num_files = length(fastq_old[[i]])
        stopifnot(num_files %% 2 == 0) # FASTQ files should come in pairs
        num_pairs = num_files / 2
        
        #   Here we take advantage of how 'list.files' orders files, and the
        #   existing naming convention for fastq_old.
        #   Specifically, mates 1 and 2 for a given pair will be adjacent, and
        #   pairs for lanes 1 and 2 will be adjacent.
        new_names = paste0(
            pd$Sample[i],
            rep(c('_R1_', '_R2_'), times = num_pairs),
            rep(1:num_pairs, each = 2),
            '.fastq.gz'
        )
        
        return(file.path(dest_dir, new_names))
    }
)

#   Create a table with columns ID, old filename, new filename. It will be in
#   "long format" in the sense that an ID will be repeated some multiple of 2
#   times. Save this table so we have a reference for how the old filenames
#   relate to the new ones

ids_long = ss(basename(unlist(fastq_new)), '_R[12]_')

file_df = data.frame(
    'sample_id' = ids_long,
    'old_path' = unlist(fastq_old),
    'new_path' = unlist(fastq_new)
)

write.csv(file_df, file = df_out_path, quote=FALSE, row.names = FALSE)

#   Now symbolically link the old FASTQs to the new destination with the new
#   naming scheme
all(file.symlink(unlist(fastq_old), unlist(fastq_new)))

session_info()
