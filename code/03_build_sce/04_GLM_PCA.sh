#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N normalize_step1_glm
#$ -o logs/normalize_step1_glm.txt
#$ -e logs/normalize_step1_glm.txt
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

## Use GPU
USAGE_CUTOFF=10
NUM_GPUS=1

avail_gpus=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader | cut -d " " -f 1 | awk -v usage="$USAGE_CUTOFF" '$1 < usage {print NR - 1}')

#  Simply exit with an error if there are no GPUs left
if [[ -z $avail_gpus ]]; then
    echo "No GPUs are available."
    exit 1
fi

export CUDA_VISIBLE_DEVICES=$(echo "$avail_gpus" | head -n $NUM_GPUS | paste -sd ",")

## Edit with your job command
Rscript normalize_step1_glm.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
