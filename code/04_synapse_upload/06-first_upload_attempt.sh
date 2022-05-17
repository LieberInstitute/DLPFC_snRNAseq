#   This script is intended to be run interactively.

#  ssh into transfer node jhpce-transfer01.jhsph.edu first, then
#  cd into the repo

#  Get absolute path to 'DLPFC_snRNAseq' repo
base_dir=$(git rev-parse --show-toplevel)

screen -S "upload_DLPFC_snRNAseq"

module load synapse/2.6.0

#   Perform a dry run as an initial test
python 06-first_upload_attempt.py \
    -c ~/synapse_credentials.yaml \
    -m "$base_dir/processed-data/04_synapse_upload/manifest.tsv" \
    -d

#   Perform the actual upload
python 06-first_upload_attempt.py \
    -c ~/synapse_credentials.yaml \
    -m "$base_dir/processed-data/04_synapse_upload/manifest.tsv" >
    $base_dir/processed-data/04_synapse_upload/06-first_upload_attempt.log 2>&1
