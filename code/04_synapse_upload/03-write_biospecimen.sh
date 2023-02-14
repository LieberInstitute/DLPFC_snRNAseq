#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G
#$ -N "write_biospecimen"
#$ -o ../../processed-data/04_synapse_upload/03-write_biospecimen.log
#$ -e ../../processed-data/04_synapse_upload/03-write_biospecimen.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
Rscript 03-write_biospecimen.R

echo "**** Job ends ****"
date
