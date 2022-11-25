library(liana)
library(tidyverse)
library(unglue)

plot_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/plots/07_ccc/"

# Per-Section lots -------------------------------------------------------
scn_fld_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/07_ccc/CCC_snRNA_R2.1/"
crn_scn <- c("Anterior", "Middle", "Posterior")
crn_scn |>
    map(.f = function(scn){
        liana_res <- readRDS(paste0(scn_fld_path, scn, "/liana_consensus.rds"))
        liana_trunc <- liana_res %>%
            # only keep interactions concordant between methods
            filter(aggregate_rank <= 0.01)

        # Create Folder
        if(!dir.exists(paste0(plot_path, scn)))
            dir.create(paste0(plot_path,scn), recursive = TRUE)

        # Save Heatmap
        jpeg(filename = paste0(plot_path, scn, "/overall.jpeg"))
        print(heat_freq(liana_trunc))
        dev.off()
    })

# Per-Sample lots -------------------------------------------------------
smp_fld_path <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/processed-data/07_ccc/CCC_snRNA_R3.1/"
file_dic <- data.frame(fl_name = list.files(smp_fld_path)) |>
    unglue_unnest(fl_name, "{Br}_{sec}", remove = FALSE) |>
    filter(substr(fl_name, start = 1, stop = 2) == "Br") |>
    left_join(
        y = data.frame(
            sec = c("ant", "mid", "post"),
            Section = c("Anterior", "Middle", "Posterior")
        )
    ) |>
    select(fl_name, Section)

map2(.x = file_dic$fl_name,
     .y = file_dic$Section,
     .f = function(name, scn){
         liana_res <- readRDS(paste0(smp_fld_path, name, "/liana_consensus.rds"))
         liana_trunc <- liana_res %>%
             # only keep interactions concordant between methods
             filter(aggregate_rank <= 0.01)

         # Create Folder
         if(!dir.exists(paste0(plot_path, scn)))
             dir.create(paste0(plot_path,scn), recursive = TRUE)

         # Save Heatmap
         jpeg(filename =  paste0(plot_path, scn, "/",name,".jpeg"))
         print(heat_freq(liana_trunc))
         dev.off()
     }
)
