
# Analysis Parameters ---------------------------------------------------
prmt <- expand.grid(
    crn_sec = c("Anterior", "Middle", "Posterior")
)

setwd("~/CCC_snRNA/")

# Helper Function for Setting Up Job --------------------------------------
start.sim <- function(
        crn_sec
) {
    #Compose job name
    job.name <- paste0("CCC_snRNA_", crn_sec)

    # NOTE:
    ## Job name has to be unique for each of your simulation settings
    ## DO NOT USE GENERIC JOB NAME FOR CONVENIENCE
    job.flag <- paste0("-N ",job.name)

    err.flag <- paste0("-e ",job.name,".err")

    out.flag <- paste0("-o ",job.name,".out")

    # Pass simulation parameters to jobs using export flag
    arg.flag <- paste0("-v crn_sec=", crn_sec)

    # Create Jobs
    system(
        paste("qsub", job.flag, err.flag, out.flag, arg.flag,
              "~/GitHub/DLPFC_snRNAseq/Boyi_R/ccc_batch_config.sh")
    )
}



# Set up job for all simulation settings ----------------------------------

for(i in 1:nrow(prmt))
    do.call(start.sim, prmt[i,,drop =FALSE])
