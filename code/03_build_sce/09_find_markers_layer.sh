#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N find_markers_layer
#$ -o logs/09_find_markers_layer.txt
#$ -e logs/09_find_markers_layer.txt
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
Rscript 09_find_markers_layer.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
