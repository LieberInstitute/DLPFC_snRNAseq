library(tidyverse)
library(unglue)


# This should be a one-time script

all_run <- c(
    "Br2720_mid", "Br2720_post", "Br2743_ant", "Br2743_mid", "Br3942_ant",
    "Br3942_mid", "Br6423_ant", "Br6423_post", "Br6432_ant", "Br6471_ant",
    "Br6471_mid", "Br6522_mid", "Br6522_post", "Br8325_ant", "Br8325_mid",
    "Br8492_mid", "Br8492_post", "Br8667_ant", "Br8667_mid"
)

runs_R3.1 <- list.files(
    path = "~/CCC_snRNA_archive/CCC_snRNA_R3.1/",
    recursive = TRUE
) |>
    unglue_data(patterns = "{Sample}/{file_name}") |>
    dplyr::filter(Sample != "log")

success_run_3.1 <- runs_R3.1 |>
    group_by(Sample) |>
    summarize(n = n()) |>
    filter(n == 2)

failed_run_R3.1 <- setdiff(all_run, success_run_3.1$Sample)
# Br8667_mid

failed_runs_R3.1 <- runs_R3.1 |>
    group_by(Sample) |>
    summarize(n = n()) |>
    filter(n != 2)


runs_R3.0 <- list.files(
    path = "~/CCC_snRNA_archive/CCC_snRNA_R3.0/",
    recursive = TRUE
) |>
    unglue_data(patterns = "{Sample}/{file_name}") |>
    dplyr::filter(Sample != "log")

success_runs_R3.0 <- runs_R3.0 |>
    group_by(Sample) |>
    summarize(n = n()) |>
    filter(n == 2)

if (length(base::intersect(success_runs_R3.0$Sample, failed_run_R3.1)) > 0) {
    map(base::intersect(success_runs_R3.0$Sample, failed_run_R3.1),
        .f = function(name) {
            folder_name <- paste0("~/CCC_snRNA_archive/CCC_snRNA_R3.1/", name, "/")
            if (!dir.exists(folder_name)) {
                dir.create(folder_name, recursive = TRUE)
            }
            file.copy(
                from = paste0("~/CCC_snRNA_archive/CCC_snRNA_R3.0/", name, "/"),
                to = paste0("~/CCC_snRNA_archive/CCC_snRNA_R3.1/"),
                recursive = TRUE
            )
        }
    )
}
