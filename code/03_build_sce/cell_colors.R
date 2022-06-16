library("here")

iwanthue_k29 <- c("#ff6785",
                 "#f50042",
                 "#a51026",
                 "#8b3929",
                 "#f57f34",
                 "#ff8400",
                 "#fdaf67",
                 "#875600",
                 "#907b00",
                 "#59532f",
                 "#bdbe64",
                 "#57a200",
                 "#a3c295",
                 "#00af70",
                 "#009e77",
                 "#008478",
                 "#0092a3",
                 "#005d79",
                 "#464a99",
                 "#bfa1ff",
                 "#9057f0",
                 "#ceade0",
                 "#822b92",
                 "#ff5ddf",
                 "#ff7bda",
                 "#ad9aa6",
                 "#ca0092",
                 "#93246a",
                 "#f30088")

iWantHue_k29 <- toupper(iwanthue_k29)

save(iWantHue_k29, file = here("processed-data", "03_build_sce","color_palletes.Rdata"))
