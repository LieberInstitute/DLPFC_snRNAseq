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


## define cell Type colors ##

cell_type_colors <- c(
  Excit = "#3264FF", #blue
  Inhib = "#D72C00", #red
  Astro = "#2C4700", # Green
  Endo.Mural = "#FF56AF", #pink
  Micro = "#4D2B70", #purple
  Micro.Oligo = "#AB0091", #magenta
  OPC = "#D2B037", # gold
  Oligo = "#E07000", #orange
  # OPC = "#AE8D00", # yellow
  drop = "black",
  Multi = "#4E586A",
  Other = "#90A583"
)

# cell_type_colors <- cell_type_colors[order(names(cell_type_colors))]
cell_type_colors <- sort(cell_type_colors)

## preview
par(las = 2) # make label text perpendicular to axis
par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
barplot(rep(1, length(cell_type_colors)),
        col = cell_type_colors,
        horiz = TRUE,
        axes = FALSE,
        names.arg = names(cell_type_colors)
)


