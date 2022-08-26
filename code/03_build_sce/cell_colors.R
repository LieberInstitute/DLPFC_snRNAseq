library("here")
library("jaffelab")

iwanthue_k29 <- c(
    "#ff6785",
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
    "#f30088"
)

iWantHue_k29 <- toupper(iwanthue_k29)

save(iWantHue_k29, file = here("processed-data", "03_build_sce", "color_palletes.Rdata"))


## define cell Type colors ##

cell_type_colors_broad <- c(
    # Excit = "#3264FF", #blue
    Excit = "#247FBC", # star command blue
    # Inhib = "#D72C00", #red
    Inhib = "#E94F37", # cinnabar
    # Astro = "#2C4700", # Green
    Oligo = "#E07000", # orange
    OPC = "#D2B037", # gold
    # OPC = "#AE8D00", # yellow
    Astro = "#3BB273", # Sea Green
    # Micro = "#4D2B70", #purple
    Micro = "#663894", # Rebecca purple
    EndoMural = "#FF56AF", # pink
    MicroOligo = "#AB0091", # magenta
    drop = "black",
    Other = "#4E586A"
)


preview_colors <- function(cell_colors) {
    par(las = 2) # make label text perpendicular to axis
    par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
    barplot(rep(1, length(cell_colors)),
        col = cell_colors,
        horiz = TRUE,
        axes = FALSE,
        names.arg = names(cell_colors)
    )
}

.scale_cell_colors <- function(color, cell_types) {
    n_ct <- length(cell_types)
    scale_colors <- grDevices::colorRampPalette(c(color, "white"))(n_ct + 1)
    scale_colors <- utils::head(scale_colors, n_ct)
    names(scale_colors) <- cell_types

    return(scale_colors)
}

expand_cell_colors <- function(cell_colors, cell_types, split = "_") {
    base_cell_types <- unique(ss(cell_types, pattern = split))
    nct <- length(base_cell_types)

    cell_colors <- cell_colors[base_cell_types]


    if (!identical(base_cell_types, cell_types)) {
        split_cell_types <- cell_types[!cell_types %in% base_cell_types]
        base_split <- rafalib::splitit(jaffelab::ss(split_cell_types, split))

        split_scale_colors <- purrr::map2(
            names(base_split), base_split,
            ~ .scale_cell_colors(
                cell_colors[[.x]],
                split_cell_types[.y]
            )
        )

        split_scale_colors <- unlist(split_scale_colors)
        cell_colors <- c(cell_colors, split_scale_colors)
    }

    return(cell_colors)
}


## All cell types
cell_types <- readLines(here("processed-data", "03_build_sce", "cell_types.txt"))
cell_type_colors <- expand_cell_colors(cell_type_colors_broad, cell_types)


## plot previews of pallets
png(here("plots", "cell_colors", "cell_colors.png"), height = 800)
preview_colors(cell_type_colors)
dev.off()


png(here("plots", "cell_colors", "cell_colors_broad.png"))
preview_colors(cell_type_colors_broad)
dev.off()

save(cell_type_colors_broad, cell_type_colors, file = here("processed-data", "03_build_sce", "cell_type_colors.Rdata"))
