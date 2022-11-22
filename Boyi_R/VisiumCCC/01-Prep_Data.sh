#!/bin/bash
## TODO: edit this section
#$ -l mem_free=10G,h_vmem=10G,h_fsize=5G
#$ -cwd
#$ -N CCC_Visium_prep
#$ -o ~/CCC_Visium/log/prep.txt
#$ -e ~/CCC_Visium/log/prep.txt
#$ -m e
#$ -M bguo6@jhu.edu

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"


## Load the R module
module load conda_R/4.2

## List current modules for reproducibility
module list

Rscript ~/GitHub/DLPFC_snRNAseq/Boyi_R/VisiumCCC/01-Prep_Data.R


echo "**** Job ends ****"
date
