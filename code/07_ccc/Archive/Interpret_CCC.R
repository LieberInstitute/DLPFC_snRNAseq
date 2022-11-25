library(liana)
library(tidyverse)



plot_path <- "/users/bguo/CCC_snRNA_archive/plots/"

scn_fld_path <- "/users/bguo/CCC_snRNA_archive/CCC_snRNA_R2.1/"
crn_scn <- c("Anterior", "Middle", "Posterior")
crn_scn |> map(
    .f = function(scn){

        liana_res <- readRDS(paste0(scn_fld_path, scn, "/liana_consensus.rds"))
        liana_trunc <- liana_res %>%
            # only keep interactions concordant between methods
            filter(aggregate_rank <= 0.01)

        # Create Folder
        if(!dir.exists(paste0(plot_path, scn)))
            dir.create(paste0(plot_path,scn), recursive = TRUE)

        # Save Heatmap
        jpeg(filename = paste0(plot_path,scn, "/overall.jpeg"))
        heat_freq(liana_trunc)
        dev.off()
    }
)

liana_test <- readRDS(paste0(scn_fld_path, "liana_test.rds"))

liana_res <- readRDS(paste0(scn_fld_path, "liana_consensus.rds"))


liana_trunc <- liana_res %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

head(liana_trunc)
head(liana_res)

heat_freq(liana_trunc)

melissa_ligand_genes <- c("SEMA6D", "FSHB", "BMP3", "EFNA5", "MDK", "DKK1", "HLA-DRA",
                          "HLA-DQA1", "HLA-DQB1", "HLA-DMA", "LRFN5", "FYN")

liana_trunc |>
    filter(ligand.complex %in% melissa_ligand_genes) |>
    chord_freq()




library(circlize)

chord_freq(liana_trunc)

liana_trunc |>
    filter(ligand.complex %in% melissa_ligand_genes) |>
    group_by(ligand.complex)

liana_res |>
    filter(receptor.complex == "FYN")

# Integrate Coronal Sections --------------------------------------------------

library(liana)
library(tidyverse)
library(here)
library(glue)

crn_sec <- c("Anterior", "Middle", "Posterior")

full_dat <- map_dfr(
    crn_sec,
    .f = function(sec){
        fl_path <- here("processed-data", "CCC", sec,
                        "liana_consensus.rds")

        readRDS(fl_path) |>
            mutate(crn_sec = sec)
    }
)

# Regardless of Source Target ---------------------------------------------

trunc_dat <- full_dat |>
    filter(aggregate_rank <= 0.01) |>
    select(crn_sec, source:receptor.complex) |>
    mutate(LR_pair = glue("{ligand.complex}-{receptor.complex}"))


# W.R.T Source Target ---------------------------------------------

trunc_src_tgt_dat <- full_dat |>
    filter(aggregate_rank <= 0.01) |>
    select(crn_sec, source:receptor.complex) |>
    mutate(
        ST_pair = glue("{source}-{target}"),
        LR_pair = glue("{ligand.complex}-{receptor.complex}"),
    )


# Most frequent LR Pair ------------------------------------------------------
trunc_dat |> group_by(LR_pair) |>
    summarize(n = n(), n_sec = n_distinct(crn_sec)) |>
    dplyr::arrange(desc(n), desc(n_sec))

# When requiring them matching ST Pair
trunc_src_tgt_dat |> group_by(LR_pair, ST_pair) |>
    summarize(n_sec = n_distinct(crn_sec)) |>
    dplyr::arrange(desc(n_sec))

# Among LR pair presented in all Sec.
trunc_src_tgt_dat |> group_by(LR_pair, ST_pair) |>
    summarize(n_sec = n_distinct(crn_sec)) |>
    dplyr::arrange(desc(n_sec)) |>
    filter(n_sec ==3) |>
    summarise(n_st_pair = n_distinct(ST_pair),
              perc_st_pair = glue("{n_st_pair}/49")) |>
    arrange(desc(n_st_pair))

# Among LR pairs that presents in all three sec
# Are they consistently prevalent in the same communication


# Venn Diagram How many are unique and how many are overlapping


# Most Popular L ------------------------------------------------------
trunc_dat |> group_by(ligand.complex) |>
    summarize(n = n(), n_sec = n_distinct(crn_sec)) |>
    dplyr::arrange(desc(n), desc(n_sec))

# Most Popular R  ------------------------------------------------------
trunc_dat |> group_by(receptor.complex) |>
    summarize(n = n(), n_sec = n_distinct(crn_sec)) |>
    dplyr::arrange(desc(n), desc(n_sec))
