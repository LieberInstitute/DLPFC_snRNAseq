#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G
#$ -N "prep_experiment"
#$ -o ../../processed-data/04_synapse_upload/01-prep_experiment.log
#$ -e ../../processed-data/04_synapse_upload/01-prep_experiment.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 01-prep_experiment.R

echo "**** Job ends ****"
date
