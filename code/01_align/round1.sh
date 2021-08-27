#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=35G,h_vmem=35G,h_fsize=100G
#$ -pe local 4
#$ -N round1
#$ -o logs/round1.$TASK_ID.txt
#$ -e logs/round1.$TASK_ID.txt
#$ -m e
#$ -t 1-3
#$ -tc 3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## load CellRanger
module load cellranger/6.1.1

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" ${JOB_NAME}.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/raw-data/FASTQ/${SAMPLE} \
    --sample=${SAMPLE} \
    --jobmode=local \
    --localcores=4 \
    --localmem=140 \
    --include-introns

## Move output
echo "Moving data to new location"
date
mkdir -p /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/cellranger/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/cellranger/

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
