#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=3G
#$ -N "write_individual"
#$ -o ../../processed-data/04_synapse_upload/04-write_individual.log
#$ -e ../../processed-data/04_synapse_upload/04-write_individual.log

echo "**** Job starts ****"
date
echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/devel
Rscript 04-write_individual.R

echo "**** Job ends ****"
date
