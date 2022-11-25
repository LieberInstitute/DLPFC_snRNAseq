
# Analysis Parameters ---------------------------------------------------
prmt <- expand.grid(
    # R2 Code
    # crn_sec = c("Anterior", "Middle", "Posterior")
    # R3 Code
    crn_sec = c("Br2720_mid", "Br2720_post", "Br2743_ant", "Br2743_mid", "Br3942_ant",
      "Br3942_mid",  "Br6423_ant",  "Br6423_post", "Br6432_ant",  "Br6471_ant",
      "Br6471_mid",  "Br6522_mid",  "Br6522_post", "Br8325_ant",  "Br8325_mid",
      "Br8492_mid",  "Br8492_post", "Br8667_ant",  "Br8667_mid" )
)


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

    err.flag <- paste0("-e /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/code/07_CCC_LIANA/logs/",job.name,".txt")

    out.flag <- paste0("-o /dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/code/07_CCC_LIANA/logs/",job.name,".txt")

    # Pass simulation parameters to jobs using export flag
    arg.flag <- paste0("-v crn_sec=", crn_sec)

    # Create Jobs
    system(
        paste("qsub", job.flag, err.flag, out.flag, arg.flag,
              "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/code/07_CCC_LIANA/01_1_ccc_batch_config.sh")
    )
}



# Set up job for all simulation settings ----------------------------------

for(i in 1:nrow(prmt))
    do.call(start.sim, prmt[i,,drop =FALSE])
