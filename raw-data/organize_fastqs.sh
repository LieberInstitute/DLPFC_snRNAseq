#!/bin/bash
#$ -cwd
#$ -t 1-15
#$ -N organize_FASTQs

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"
echo "****"

# Iteratively create soft links for FASTQ files, referencing a TSV with library sample name/some actual sample info

SAMPLEID=$(awk 'BEGIN {FS="\t"} {print $1}' ./sample_libs_rounds2-5.tsv | awk "NR==${SGE_TASK_ID}")

mkdir FASTQ/${SAMPLEID}/
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-11-22_KMay110521/${SAMPLEID}_L00*/* FASTQ/${SAMPLEID}/

echo "**** Job completed ****"
date
