#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -pe local 4
#$ -N round0
#$ -o logs/round0.txt
#$ -e logs/round0.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"


## List current modules for reproducibility
module list

## Load Cell Ranger
module load cellranger/6.1.1

## Locate file
SAMPLE=$(awk "NR==1" ${JOB_NAME}.txt)
echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A \
    --fastqs=/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/raw-data/FASTQ/${SAMPLE} \
    --sample=Br2743-DLPFC-mid \
    --jobmode=local \
    --localcores=4 \
    --localmem=20 \
    --include-introns

## Move output
echo "Moving data to new location"
date
#mkdir -p /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/cellranger/
mv ${SAMPLE} /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/cellranger/


echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/
