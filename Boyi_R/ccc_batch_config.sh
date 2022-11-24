#!/bin/bash
## TODO: edit this section
#$ -l mem_free=45G,h_vmem=45G,h_fsize=5G
#$ -cwd
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


# Create a directory to store your log files
# Import to debuging
resPath=~/CCC_snRNA
logPath=~/CCC_snRNA/log/
mkdir -p $logPath $resPath

# Run simulation with simulation parameters using R batch --args flag
###      The $ sign helps to fetch the simulation parameters passed from start_sim.R\
###      If the parameter is a string (a vector of string), you need to quote the string with single quote `
R CMD BATCH --no-save "--args crn_sec='$crn_sec'" ~/GitHub/DLPFC_snRNAseq/Boyi_R/CCC_with_LIANA.R  $logPath${crn_sec}.Rout


echo "**** Job ends ****"
date
