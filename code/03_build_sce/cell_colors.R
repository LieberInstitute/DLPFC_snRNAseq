library("here")
library("jaffelab")
library("RColorBrewer")
library("SingleCellExperiment")

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
    Ambiguous = "#A9998F", # light brown-grey
    drop = "black"
    # Other = "#4E586A",
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
(cell_type_colors <- cell_type_colors[cell_types])
# Astro    EndoMural EndoMural_01 EndoMural_02        Micro   MicroOligo     Oligo_01     Oligo_02     Oligo_03
# "#3BB273"    "#FF56AF"    "#FF56AF"    "#FFAAD7"    "#663894"    "#AB0091"    "#E07000"    "#EA9F55"    "#F4CFAA"
# OPC     Excit_01     Excit_02     Excit_03     Excit_04     Excit_05     Excit_06     Excit_07     Excit_08
# "#D2B037"    "#247FBC"    "#3287C0"    "#4190C4"    "#4F98C9"    "#5EA1CD"    "#6DA9D2"    "#7BB2D6"    "#8ABADB"
# Excit_09     Excit_10     Excit_11     Excit_12     Excit_13     Excit_14     Excit_15     Inhib_01     Inhib_02
# "#98C3DF"    "#A7CBE4"    "#B5D4E8"    "#C4DCED"    "#D3E5F1"    "#E1EDF6"    "#F0F6FA"    "#E94F37"    "#EC6C58"
# Inhib_03     Inhib_04     Inhib_05     Inhib_06    Ambiguous      drop_01      drop_02      drop_03      drop_04
# "#F08979"    "#F3A79B"    "#F7C4BC"    "#FBE1DD"    "#A9998F"    "#000000"    "#3F3F3F"    "#7F7F7F"    "#BFBFBF"

## plot previews of pallets
png(here("plots", "cell_colors", "cell_colors.png"), height = 800)
preview_colors(cell_type_colors)
dev.off()


png(here("plots", "cell_colors", "cell_colors_broad.png"))
preview_colors(cell_type_colors_broad)
dev.off()

save(cell_type_colors_broad, cell_type_colors, file = here("processed-data", "03_build_sce", "cell_type_colors.Rdata"))

## Add colors for layer annotation cell types
library("SingleCellExperiment")
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)
(cell_types_layer_all <- levels(sce$cellType_layer))
# [1] "Astro"        "EndoMural"    "Micro"        "Oligo"        "OPC"          "Excit_L2/3"   "Excit_L3"
# [8] "Excit_L3/4/5" "Excit_L4"     "Excit_L5"     "Excit_L5/6"   "Excit_L6"     "Excit_ambig"  "Inhib"

cell_types_layer <- cell_types_layer_all[cell_types_layer_all != "Excit_ambig"]

cell_type_colors_layer <- expand_cell_colors(cell_type_colors_broad, cell_types_layer)

png(here("plots", "cell_colors", "cell_type_colors_layer.png"))
preview_colors(cell_type_colors_layer)
dev.off()


## hmm blues are very simmilar...try brewer blues
blues <- brewer.pal(9, "PuBu")

cell_type_colors_layer[grepl("Excit_", names(cell_type_colors_layer))] <- blues[3:9]

cell_type_colors_layer <- cell_type_colors_layer[cell_types_layer]

png(here("plots", "cell_colors", "cell_type_colors_layer_brew.png"))
preview_colors(cell_type_colors_layer)
dev.off()

## In to green gradiant
# blues <- brewer.pal(8, "YlGnBu")
# cell_type_colors_layer[grepl("Excit_", names(cell_type_colors_layer))] <- blues[3:8]
#
# png(here("plots", "cell_colors", "cell_type_colors_layer_brew.png"))
# preview_colors(cell_type_colors_layer)
# dev.off()

cell_type_colors_layer <- c(cell_type_colors_layer, Excit_ambig = "#A0A7A7")
cell_type_colors_layer <- cell_type_colors_layer[cell_types_layer_all]
# Astro    EndoMural        Micro        Oligo          OPC   Excit_L2/3     Excit_L3 Excit_L3/4/5     Excit_L4
# "#3BB273"    "#FF56AF"    "#663894"    "#E07000"    "#D2B037"    "#D0D1E6"    "#A6BDDB"    "#74A9CF"    "#3690C0"
# Excit_L5   Excit_L5/6     Excit_L6  Excit_ambig        Inhib
# "#0570B0"    "#045A8D"    "#023858"    "#A0A7A7"    "#E94F37"

## Save
save(cell_type_colors_layer, file = here("processed-data", "03_build_sce", "cell_type_colors_layer.Rdata"))

## Save in metadata
metadata(sce)$cell_type_colors_layer <- cell_type_colors_layer
save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))
