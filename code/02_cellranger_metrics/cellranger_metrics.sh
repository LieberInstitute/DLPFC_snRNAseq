#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N cellranger_metrics
#$ -o logs/cellranger_metrics.txt
#$ -e logs/cellranger_metrics.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 02_cellranger_metrics.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
