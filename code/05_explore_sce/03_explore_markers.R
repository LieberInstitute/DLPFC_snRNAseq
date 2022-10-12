library("SingleCellExperiment")
# library("scater")
library("tidyverse")
library("EnhancedVolcano")
library("ComplexHeatmap")
library("here")
library("sessioninfo")

#### Load data ####
load(here("processed-data", "sce", "sce_DLPFC.Rdata"), verbose = TRUE)

## markers data
load(here("processed-data", "03_build_sce", "cell_type_markers.Rdata"), verbose = TRUE)
# markers_1vALL ## findMarkers - all cell types
# markers_mean_ratio ## Mean ratio - all cell types
# markers_mean_ratio_broad

load(here("processed-data", "03_build_sce", "cell_type_markers_1vALL_mod.Rdata"), verbose = TRUE)
# markers_1vALL_sample ## enrichment - all cell types
# Missing enrichment - cell subtype - TODO

load(here("processed-data", "03_build_sce", "subtype_markers.Rdata"), verbose = TRUE)
# subtype_markers ## findMarkers - cell subtype - 'All' and "Any' settings

#### set up plotting ####
source(here("code", "03_build_sce", "my_plotExpression.R"))

plot_dir <- here("plots", "05_explore_sce", "03_explore_markers")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# cell_type_colors_broad
cell_type_colors <- metadata(sce)$cell_type_colors[levels(sce$cellType_hc)]

#### 1vALL markers ####
sapply(markers_1vALL, function(x) sum(x$FDR < 0.05))
# Astro Endo.Mural_01 Endo.Mural_02      Excit_01      Excit_02      Excit_03      Excit_04      Excit_05
# 482           914            13            53            21            81            42             0
# Excit_06      Excit_07      Excit_08      Excit_09      Excit_10      Excit_11      Excit_12      Excit_13
# 237           185            86             0           742           103           144           442
# Excit_14      Excit_15      Inhib_01      Inhib_02      Inhib_03      Inhib_04      Inhib_05      Inhib_06
# 0             0           240           131             5           178            97            74
# Micro      Oligo_01      Oligo_02      Oligo_03           OPC
# 1016            13           889             0           277

(markers_1vALL_top10 <- sapply(markers_1vALL, function(x) rownames(x[1:10, ])))
#       Astro        Endo.Mural_01 Endo.Mural_02 Excit_01     Excit_02     Excit_03     Excit_04     Excit_05
# [1,] "AC012405.1" "ATP10A"      "COL6A3"      "AC095050.1" "AC011754.2" "AL450352.1" "CNGB3"      "AL121929.1"
# [2,] "ITGB4"      "ABCB1"       "TAGLN"       "AL157944.1" "LINC01435"  "COL22A1"    "AP003066.1" "AC026474.1"
# [3,] "SLC16A9"    "FLT1"        "BMP5"        "AC009975.2" "PDE4D"      "CCDC168"    "CASC15"     "AC090371.2"
# [4,] "RNF182"     "COBLL1"      "MYH11"       "CUX2"       "SORCS1"     "AL356295.1" "LINC01680"  "ANXA8L1"
# [5,] "TNC"        "LEF1"        "KCNK17"      "LINC02015"  "AC099517.1" "DCC"        "AC073091.3" "AC061975.4"
# [6,] "AL136366.1" "THSD4"       "COL1A1"      "AC117453.1" "HS3ST2"     "LINC01725"  "AC016687.2" "AC111182.2"
# [7,] "BBOX1"      "SLC7A1"      "ADAMTS1"     "AC016042.1" "AC007368.1" "AL353595.1" "AL158835.1" "TRAV38-1"
# [8,] "EDNRB"      "RGS5"        "SRPX2"       "LINC01378"  "AC011754.1" "AC073365.1" "ATP12A"     "AC007278.1"
# [9,] "AL627316.1" "GPCPD1"      "ITK"         "PLEKHS1"    "ARHGEF28"   "SLC22A10"   "TRABD2A"    "LINC01981"
# [10,] "AC103874.1" "EPAS1"       "AC099782.2"  "NDST3"      "IL17REL"    "ABCA4"      "AC016687.3" "IFNA4"
#       Excit_06     Excit_07     Excit_08     Excit_09     Excit_10 Excit_11     Excit_12     Excit_13     Excit_14
# [1,] "SMYD1"      "NPSR1-AS1"  "LINC02718"  "IDS"        "YWHAB"  "NTNG1"      "TRPM3"      "AC132872.1" "AC129926.2"
# [2,] "ZNF804B"    "FER1L6-AS2" "LINC02726"  "FBXL16"     "CALM1"  "FOXP2"      "GPR151"     "MIR4458HG"  "SLC12A8"
# [3,] "AL136119.1" "ITGA8"      "AL358335.2" "ADRA2A"     "STMN2"  "AC019068.1" "KLHL1"      "MYLPF"      "CCDC196"
# [4,] "OLFML2B"    "IFNG-AS1"   "AC002331.1" "INSYN1"     "CREG2"  "SCN7A"      "POU4F1"     "PHPT1"      "AC023824.3"
# [5,] "SORBS2"     "GRM4"       "NPFFR2"     "MTURN"      "MAP1B"  "SULF1"      "TMEM163"    "ATP5ME"     "BX890604.2"
# [6,] "NTNG2"      "AC068587.4" "FRMPD4"     "DDN"        "LMO4"   "CYP27C1"    "AC233296.1" "GFAP"       "AC124290.1"
# [7,] "MCTP2"      "TRMT9B"     "AC104689.2" "SOWAHA"     "BEX1"   "COL6A5"     "RASGEF1B"   "MIF-AS1"    "TSIX"
# [8,] "LINC02197"  "ASIC2"      "HS3ST4"     "AP002449.1" "KCNA2"  "SLITRK6"    "CDH8"       "SAT2"       "CNTD1"
# [9,] "RGS12"      "LINC01885"  "MCTP1"      "RPL22"      "SCN8A"  "EMB"        "PRKD1"      "ARL2"       "CACNG3"
# [10,] "KCTD16"     "GHR"        "LINC01606"  "FP565260.6" "BASP1"  "TCF7L2"     "LINC01497"  "MIF"        "AC005972.3"
#       Excit_15     Inhib_01     Inhib_02    Inhib_03     Inhib_04     Inhib_05     Inhib_06     Micro
# [1,] "GEM"        "OTX2-AS1"   "DPP10-AS3" "TFAP2D"     "AC137770.1" "AL391832.4" "NMU"        "APBB1IP"
# [2,] "HSPB1"      "OTX2"       "MYO5B"     "EBF3"       "AC132803.1" "PROX1"      "FLT3"       "ADAM28"
# [3,] "AC016722.2" "LINC01210"  "ZNF804A"   "NR5A2"      "PRELID2"    "PPP1R1C"    "AC113347.4" "AC008691.1"
# [4,] "HSPD1"      "AL161757.4" "CNTNAP3B"  "VWA5B1"     "LINC01344"  "PROX1-AS1"  "CDH9"       "FYB1"
# [5,] "GTPBP1"     "GATA3"      "LRRC38"    "SIM1"       "ROR2"       "KMO"        "AC073332.1" "DLEU1"
# [6,] "AC008105.1" "CASR"       "ERBB4"     "AC010086.3" "PKP2"       "CXCL14"     "AC011586.2" "CSF2RA"
# [7,] "HSPA9"      "TSPAN18"    "ANK1"      "LHX1-DT"    "PTCHD4"     "VIP"        "GRIK3"      "LINC01374"
# [8,] "WNT10B"     "AC007159.1" "SCN1A-AS1" "KRT8"       "HAPLN1"     "CHRNA2"     "AC076968.2" "SYK"
# [9,] "YARS"       "SOX14"      "NHS"       "NMD3"       "PDGFD"      "AL450338.2" "TRIM67"     "CX3CR1"
# [10,] "SRSF3"      "AC245187.2" "CNTNAP3"   "PITX2"      "FGF13"      "AC020930.1" "SATB1-AS1"  "LNCAROD"
#       Oligo_01    Oligo_02      Oligo_03     OPC
# [1,] "MTRNR2L12" "SYNJ2"       "PAIP2B"     "AC004852.2"
# [2,] "MT-ND3"    "SHROOM4"     "AC026401.2" "SMOC1"
# [3,] "MT-CO3"    "AC008571.2"  "NKAP"       "FERMT1"
# [4,] "MT-CO1"    "FAM107B"     "AC090515.5" "STK32A"
# [5,] "MT-ND2"    "PRR5L"       "TMEM95"     "COL9A1"
# [6,] "MT-CO2"    "MOG"         "LINC00463"  "PDGFRA"
# [7,] "MT-CYB"    "LINC00639"   "FCRLA"      "AL512308.1"
# [8,] "MT-ND5"    "SLC7A14-AS1" "FAM225A"    "BEST3"
# [9,] "MT-ATP6"   "LRP2"        "SLC26A9"    "AL445250.1"
# [10,] "MT-ND4"    "LINC01170"   "AC090503.1" "CACNG4"

sapply(markers_1vALL, function(x) sum(x$FDR < 0.05) == 0)

#### Plot 1vALL markers ####
markers_1vALL_top10 <- as.list(as.data.frame(markers_1vALL_top10))

my_plotMarkers(
    sce = sce,
    marker_list = markers_1vALL_top10,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_1vALL_top10.pdf")
)

## volcano plots
pdf(here(plot_dir, "markers_1vALL_volcano.pdf"), width = 10, height = 10)
for (n in names(markers_1vALL)) {
    message(n)

    markers <- markers_1vALL[[n]]
    st <- paste0("n nuc=", sum(sce$cellType_hc == n), ", n genes FDR<0.05=", sum(markers$FDR < 0.05))
    message(st)
    print(EnhancedVolcano(markers,
        lab = rownames(markers),
        x = "summary.logFC",
        y = "FDR",
        pCutoff = 0.05,
        xlab = bquote(~ Log[2] ~ "summary fold change"),
        ylab = bquote(~ -Log[10] ~ italic(FDR)),
        title = n,
        subtitle = st
    ))
}
dev.off()

#### Bind in to one big table ####
markers_1vALL_df <- do.call(
    "rbind",
    lapply(markers_1vALL, function(x) as.data.frame(x[, c("p.value", "FDR", "summary.logFC")]))
) |>
    rename_with(.fn = ~ paste0(.x, "_1vALL")) |>
    rownames_to_column("ct_gene") |>
    mutate(ct_gene = gsub("Endo.Mural", "EndoMural", ct_gene)) |>
    separate(ct_gene, into = c("cellType.target", "gene"), extra = "merge", sep = "\\.")

head(markers_1vALL_df)
#   cellType.target       gene p.value_1vALL     FDR_1vALL summary.logFC_1vALL
# 1           Astro AC012405.1 2.229583e-118 8.160497e-114           0.4273925
# 2           Astro      ITGB4 2.015016e-110 3.687580e-106           0.2189336
# 3           Astro    SLC16A9  1.148270e-81  1.400928e-77           0.2557755
# 4           Astro     RNF182  1.981230e-81  1.812875e-77           0.1530592
# 5           Astro        TNC  1.329933e-61  9.735377e-58           0.2517694
# 6           Astro AL136366.1  3.721927e-61  2.270438e-57           0.1070652

## should have 36601 gens for each cell type
markers_1vALL_df |> count(cellType.target)

markers_1vALL_df |>
    group_by(cellType.target) |>
    summarize(n_markers = sum(FDR_1vALL < 0.05)) |>
    arrange(n_markers)

#   cellType.target n_markers
# 1 Excit_05                0
# 2 Excit_09                0
# 3 Excit_14                0
# 4 Excit_15                0
# 5 Oligo_03                0

#### Mean Ratio Markers ####
markers_mean_ratio_top10 <- markers_mean_ratio |>
    group_by(cellType.target) |>
    slice(1:10) |>
    select(gene, cellType.target) |>
    unstack() |>
    as.list()

## Cell types with out good markers
markers_mean_ratio |>
    group_by(cellType.target) |>
    summarize(n_markers = sum(ratio > 1)) |>
    arrange(n_markers)

# cellType.target n_markers
# <fct>               <int>
# 1 Excit_05                0 *
# 2 Oligo_03                0 *
# 3 Endo.Mural_02           1
# 4 Excit_15                1 *
# 5 Excit_09                4 *
# 6 Inhib_03                6

markers_mean_ratio |> filter(cellType.target == "Excit_09")

my_plotMarkers(
    sce = sce,
    marker_list = markers_mean_ratio_top10,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_mean_ratio_top10.pdf")
)


markers_mean_ratio_broad_top10 <- markers_mean_ratio_broad |>
    group_by(cellType.target) |>
    slice(1:10) |>
    select(gene, cellType.target) |>
    unstack() |>
    as.list()

my_plotMarkers(
    sce = sce,
    marker_list = markers_mean_ratio_broad_top10,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_mean_ratio_broad_top10.pdf")
)

my_plotMarkers(
    sce = sce,
    marker_list = markers_mean_ratio_broad_top10,
    cat = "cellType_broad_hc",
    fill_colors = cell_type_colors_broad,
    pdf_fn = here(plot_dir, "markers_mean_ratio_broad_broad_top10.pdf")
)


markers_mean_ratio |>
    filter(gene == "MALAT1", rank_ratio <= 25) |>
    select(gene, cellType.target, ratio, rank_ratio)

# gene   cellType.target ratio rank_ratio
# <chr>  <fct>           <dbl>      <int>
# 1 MALAT1 Excit_05        0.837          4 * poor markers for 1vALL and MR
# 2 MALAT1 Endo.Mural_02   0.965          2
# 3 MALAT1 Oligo_03        0.846         14 *
# 4 MALAT1 Excit_15        0.829          4

## hm - problem cell types?
markers_mean_ratio |>
    filter(ratio < 1, rank_ratio <= 25) |>
    count(cellType.target)
# # A tibble: 7 × 2
# cellType.target     n
# <fct>           <int>
#   1 Endo.Mural_02      24
# 2 Excit_05           25
# 3 Excit_09           21
# 4 Excit_15           24
# 5 Inhib_03           19
# 6 Oligo_01           12
# 7 Oligo_03           25

head(markers_1vALL_df)

markers_all <- markers_mean_ratio |>
    left_join(markers_1vALL_df) |>
    mutate(cellType.target = factor(cellType.target))

markers_all |> count(cellType.target)

## hmmm funky because of summary.logFC?
ratio_vs_stdFC <- markers_all |>
    ggplot(aes(x = ratio, y = summary.logFC_1vALL)) +
    geom_point() +
    facet_wrap(~cellType.target, scales = "free_x") +
    theme_bw()

ggsave(ratio_vs_stdFC, filename = here(plot_dir, "ratio_vs_stdFC.png"), height = 12, width = 12)

## volcano plots
pdf(here(plot_dir, "markers_mean_ratio_volcano.pdf"), width = 10, height = 10)
for (n in levels(markers_all$cellType.target)[1:3]) {
    markers <- markers_all |> filter(cellType.target == n)
    message(n, nrow(markers))
    st <- paste0("n nuc=", sum(sce$cellType_hc == n), ", n genes FDR<0.05=", sum(markers$FDR_1vALL < 0.05))
    message(st)

    # print(EnhancedVolcano(markers,
    #                       lab = rownames(markers),
    #                       x = 'ratio',
    #                       y = 'FDR_1vALL',
    #                       pCutoff = 0.05,
    #                       xlab = "Mean Ratio",
    #                       ylab = bquote(~-Log[10] ~ italic(FDR)),
    #                       title = n,
    #                       subtitle = st))
}
dev.off()

#### Explore 1vALL control for Sample ####
markers_1vALL_sample |>
    filter(log.FDR < log(0.05)) |>
    count(cellType.target) |>
    arrange(n)

pdf(here(plot_dir, "Volcano_markers_1vALL_Sample.pdf"), width = 10, height = 10)
for (n in levels(markers_1vALL_sample$cellType.target)) {
    message(n)

    markers <- markers_1vALL_sample |>
        filter(cellType.target == n) |>
        mutate(FDR = exp(log.FDR))

    st <- paste0("n nuc=", sum(sce$cellType_hc == n), ", n genes FDR<0.05=", sum(markers$FDR < 0.05))

    message(st)
    print(EnhancedVolcano(markers,
        lab = markers$gene,
        x = "logFC",
        y = "FDR",
        pCutoff = 0.05,
        xlab = bquote(~ Log[2] ~ "summary fold change"),
        ylab = bquote(~ -Log[10] ~ italic(FDR)),
        title = n,
        subtitle = st
    ))
}
dev.off()

## How does custom 1vALL compare to standard 1vALL? ###

## cell types with poor markers
markers_1vALL_sample |>
    group_by(cellType.target) |>
    summarize(n_markers = sum(log.FDR < log10(0.05))) |>
    arrange(n_markers)

# cellType.target n_markers ## looks better than standard 1vALL?
# <fct>               <int>
#   1 Excit_14             1058
# 2 Excit_15             1348
# 3 Excit_11             1973
# 4 Excit_12             2865
# 5 Oligo_01             2928

markers_1vALL_sample |>
    group_by(cellType.target) |>
    summarize(n_markers = sum(logFC > 0)) |>
    arrange(n_markers)


markers_1vALL_sample_top10 <- markers_1vALL_sample |>
    group_by(cellType.target) |>
    slice(1:10) |>
    select(gene, cellType.target) |>
    unstack() |>
    as.list()

map2(markers_1vALL_sample_top10, markers_1vALL_top10, intersect)


my_plotMarkers(
    sce = sce,
    marker_list = markers_1vALL_sample_top10,
    cat = "cellType_hc",
    fill_colors = cell_type_colors,
    pdf_fn = here(plot_dir, "markers_mean_ratio_top10.pdf")
)


compare_1vALL <- markers_1vALL_sample |>
    select(gene, cellType.target, logFC_Sample = logFC, std.logFC_Sample = std.logFC, log.FDR_Sample = log.FDR) |>
    mutate(FDR_Sample = exp(log.FDR_Sample)) |>
    inner_join(markers_1vALL_df) |>
    mutate(signif = case_when(
        FDR_1vALL < 0.05 & log.FDR_Sample < log(0.05) ~ "Signif_both",
        FDR_1vALL < 0.05 ~ "Signif_1vALL",
        log.FDR_Sample < log(0.05) ~ "Signif_Sample",
        TRUE ~ "None"
    ))

compare_1vALL |> count(signif)

scatter_logFC_1vALL <- compare_1vALL |>
    ggplot(aes(x = summary.logFC_1vALL, y = logFC_Sample, color = signif)) +
    geom_point(size = 1, alpha = .3) +
    facet_wrap(~cellType.target)

ggsave(scatter_logFC_1vALL, filename = here(plot_dir, "scatter_logFC_1vALL.png"), height = 12, width = 12)


scatter_FDR_1vALL <- compare_1vALL |>
    ggplot(aes(x = log(FDR_1vALL), y = log.FDR_Sample, color = signif)) +
    geom_point(size = 1, alpha = .3) +
    facet_wrap(~cellType.target)

ggsave(scatter_FDR_1vALL, filename = here(plot_dir, "scatter_FDR_1vALL.png"), height = 12, width = 12)

#### Compare to mean ratio ####
markers_all_sample <- markers_mean_ratio |>
    left_join(markers_1vALL_sample)

markers_all_sample |> count(cellType.target)

## hmmm funky because of summary.logFC?
ratio_vs_stdFC_sample <- markers_all_sample |>
    ggplot(aes(x = ratio, y = std.logFC)) +
    geom_point(size = 0.5, alpha = 0.5) +
    facet_wrap(~cellType.target) +
    theme_bw() +
    xlim(0, 25) ## cuts off nice points, but need to zoom

ggsave(ratio_vs_stdFC_sample, filename = here(plot_dir, "ratio_vs_stdFC_Sample.png"), height = 12, width = 12)


#### subtype Markers ####
sapply(markers_1vALL, function(x) sum(x$FDR < 0.05))

names(subtype_markers$Excit)
names(subtype_markers)

subtype_markers <- transpose(subtype_markers)

n_markers_all <- map(subtype_markers$all, ~ sapply(.x, function(x) sum(x$FDR < 0.05)))
map(n_markers_all, ~ .x[!is.na(.x)])
# $EndoMural
# EndoMural_01 EndoMural_02
# 3107         1853
#
# $Excit
# Excit_01 Excit_02 Excit_03 Excit_04 Excit_05 Excit_06 Excit_07 Excit_08 Excit_09 Excit_10 Excit_11 Excit_12 Excit_13
# 79       35      144       65        2      391      368      179        1      834      264      421      843
# Excit_14 Excit_15
# 0        0
#
# $Inhib
# Inhib_01 Inhib_02 Inhib_03 Inhib_04 Inhib_05 Inhib_06
# 1258     1429      413      909      774      491
#
# $Oligo
# Oligo_01 Oligo_02 Oligo_03
# 9268     5296       83

n_markers_any <- map(subtype_markers$any, ~ sapply(.x, function(x) sum(x$FDR < 0.05)))
# map(n_markers_any,~.x[!is.na(.x)])
# $EndoMural
# EndoMural_01 EndoMural_02
# 3107         1853
#
# $Excit
# Excit_01 Excit_02 Excit_03 Excit_04 Excit_05 Excit_06 Excit_07 Excit_08 Excit_09 Excit_10 Excit_11 Excit_12 Excit_13
# 29758    25857    21597    26083    18577    13982    14395    22227    11218     9565    12606    12215    15609
# Excit_14 Excit_15
# 3559      270
#
# $Inhib
# Inhib_01 Inhib_02 Inhib_03 Inhib_04 Inhib_05 Inhib_06
# 6772    14122    11326    10240    13736    12297
#
# $Oligo
# Oligo_01 Oligo_02 Oligo_03
# 13169     8438     9270


subtype_markers$all$Excit

subtype_markers$all$Excit$Excit_04
subtype_markers$all$Excit$Excit_05
subtype_markers$all$Excit$Excit_09
subtype_markers$all$Excit$Excit_14

sapply(subtype_markers$all$Excit, function(x) rownames(x[is.na(x$FDR)])[[1]])


sum(subtype_markers$all$Excit$Excit_15$FDR < 0.05)


#### Heatmaps ####

mean_ratio_top <- markers_mean_ratio |>
    group_by(cellType.target) |>
    slice(1:2) |>
    filter(gene != "MALAT1") |>
    mutate(cellType.target = gsub("\\.", "", cellType.target))

mean_ratio_top |>
    ungroup() |>
    count(gene) |>
    filter(n > 1)

## Pseudobulk by cell type
pb_sce_hc <- scuttle::aggregateAcrossCells(sce,
    ids = sce$cellType_hc
)

## Need to normalize?
sizeFactors.PB <- scater::librarySizeFactors(assays(pb_sce_hc)$counts)

# # Normalize with these LSFs
assays(pb_sce_hc)$logcounts <- t(apply(assays(pb_sce_hc)$counts, 1, function(x) {
    log2(x / sizeFactors.PB + 1)
}))


summary(colSums(assays(pb_sce_hc)$logcounts))


marker_anno <- mean_ratio_top |>
    column_to_rownames("gene") |>
    select(cellType = cellType.target)

ct_anno <-
    png(here(plot_dir, "heatmap_test_pheatmap.png"), height = 800, width = 800)
# pheatmap(assays(pb_sce_hc_top)$counts,
pheatmap(assays(pb_sce_hc)$logcounts[mean_ratio_top$gene, ],
    main = "Top Markers",
    annotation_row = marker_anno,
    annotation_colors = list(cellType = cell_type_colors)
)
dev.off()

png(here(plot_dir, "heatmap_test_complex.png"), height = 800, width = 800)
Heatmap(assays(pb_sce_hc_top)$counts)
dev.off()


#### Layer markers ####
## explore markers
load(here("processed-data", "03_build_sce", "cell_type_markers_layer.Rdata"), verbose = TRUE)
# markers_1vALL
# markers_1vALL_enrich
# markers_mean_ratio

markers_1vALL_enrich |>
    filter(log.FDR < log(0.05)) |>
    group_by(cellType.target) |>
    count()

# cellType.target     n
# <fct>           <int>
#   1 Astro            7083
# 2 EndoMural        7981
# 3 Micro            6616
# 4 Oligo            4258
# 5 OPC              5107
# 6 Excit_L2/3       1063
# 7 Excit_L3        10114
# 8 Excit_L3/4/5     6226
# 9 Excit_L4         7300
# 10 Excit_L5         7645
# 11 Excit_L5/6       6872
# 12 Excit_L6         6317
# 13 Inhib            9019


markers_mean_ratio |>
    filter(ratio > 1) |>
    group_by(cellType.target) |>
    count()

# cellType.target     n
# <fct>           <int>
# 1 Astro             242
# 2 EndoMural         261
# 3 Micro             305
# 4 Oligo              27
# 5 OPC               353
# 6 Excit_L2/3        184
# 7 Excit_L3          190
# 8 Excit_L3/4/5      290
# 9 Excit_L4          288
# 10 Excit_L5          364
# 11 Excit_L5/6        425
# 12 Excit_L6          667
# 13 Inhib             357

layer_markers <- left_join(markers_1vALL_enrich, markers_mean_ratio)

## version from spatialDLPFC from Nick
marker_stats_spatial_layer <- readRDS(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/05-shared_utilities/marker_stats_layer.rds")
marker_stats_spatial_broad <- readRDS(file = "/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/05-shared_utilities/marker_stats_broad.rds")
head(marker_stats_spatial)

## do stats match?
all_equal(markers_mean_ratio[colnames(markers_mean_ratio)], markers_mean_ratio)
# [1] TRUE ## nice get_mean_ratio2 is consistent :)

## Nicks final list
marker_list_layer <- readLines("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/05-shared_utilities/markers_layer.txt")
marker_list_broad <- readLines("/dcs04/lieber/lcolladotor/spatialDLPFC_LIBD4035/spatialDLPFC/processed-data/spot_deconvo/05-shared_utilities/markers_broad.txt")

## layer markers
layer_marker_stats_anno <- marker_stats_spatial_layer |>
    add_column(symbol = rownames(sce)[match(marker_stats_spatial_layer$gene, rowData(sce)$gene_id)]) |>
    filter(!grepl("^MT-", symbol)) |> ## Exclude mitochondrial genes
    mutate(marker = gene %in% marker_list_layer) |>
    filter(rank_ratio <= 50 & marker) |> ## this catches the 20 Oligo markers, not all top 25
    select(symbol, gene, cellType.target, rank_ratio) |>
    group_by(cellType.target) |>
    mutate(adj_rank_ratio = row_number()) |> ## adjust for removed MT genes
    mutate(marker_anno = paste0(cellType.target, "-", adj_rank_ratio))


setequal(layer_marker_stats_anno$gene, marker_list_layer)
# [1] TRUE
layer_marker_stats_anno |> count(cellType.target)
layer_marker_stats_anno |>
    group_by(cellType.target) |>
    count(rank_ratio == adj_rank_ratio)

## Add to rowdata
rowData(sce)$marker_layer <- layer_marker_stats_anno$marker_anno[match(rowData(sce)$gene_id, layer_marker_stats_anno$gene)]


## braod markers
broad_marker_stats_anno <- marker_stats_spatial_broad |>
    add_column(symbol = rownames(sce)[match(marker_stats_spatial_broad$gene, rowData(sce)$gene_id)]) |>
    filter(!grepl("^MT-", symbol)) |>
    mutate(marker = gene %in% marker_list_broad) |>
    filter(rank_ratio <= 50 & marker) |> ## this catches the 20 Oligo markers, not all top 25
    select(symbol, gene, cellType.target, rank_ratio) |>
    group_by(cellType.target) |>
    mutate(adj_rank_ratio = row_number()) |> ## adjust for removed MT genes
    mutate(marker_anno = paste0(cellType.target, "-", adj_rank_ratio))

setequal(broad_marker_stats_anno$gene, marker_list_broad)
# [1] TRUE
broad_marker_stats_anno |> count(cellType.target)
broad_marker_stats_anno |>
    group_by(cellType.target) |>
    count(rank_ratio == adj_rank_ratio)

rowData(sce)$marker_broad_spatial <- broad_marker_stats_anno$marker_anno[match(rowData(sce)$gene_id, broad_marker_stats_anno$gene)]

layer_marker_stats_anno |> filter(rank_ratio == 1)

## check some genes
rowData(sce)[c("GRIP2", "MBP", "SLC25A18", "VCAN"), ]
# DataFrame with 4 rows and 9 columns
# source     type         gene_id gene_version   gene_name      gene_type binomial_deviance marker_layer
# <factor> <factor>     <character>  <character> <character>    <character>         <numeric>  <character>
#   GRIP2      HAVANA     gene ENSG00000144596           13       GRIP2 protein_coding           84376.8     Inhib-12
# MBP        HAVANA     gene ENSG00000197971           15         MBP protein_coding          962309.0      Oligo-3
# SLC25A18   HAVANA     gene ENSG00000182902           14    SLC25A18 protein_coding           88904.8     Astro-11
# VCAN       HAVANA     gene ENSG00000038427           16        VCAN protein_coding          232054.0       OPC-17
# marker_broad_spatial
# <character>
#   GRIP2                 Inhib-2
# MBP                   Oligo-6
# SLC25A18             Astro-12
# VCAN                   OPC-23

## Save annotated rowData
# save(sce, file = here("processed-data", "sce", "sce_DLPFC.Rdata"))


#### Markere Gene Heatmap ####
pb_sce_layer <- scuttle::aggregateAcrossCells(sce[, !is.na(sce$cellType_layer)],
    ids = sce$cellType_layer[!is.na(sce$cellType_layer)]
)

## Need to normalize?
sizeFactors.PB <- scater::librarySizeFactors(assays(pb_sce_layer)$counts)

# # Normalize with these LSFs
assays(pb_sce_layer)$logcounts <- t(apply(assays(pb_sce_layer)$counts, 1, function(x) {
    log2(x / sizeFactors.PB + 1)
}))

layer_markers_top <- markers_mean_ratio |>
    filter(rank_ratio <= 2) |>
    arrange(cellType.target)

marker_anno <- layer_markers_top |>
    column_to_rownames("gene") |>
    select(cellType = cellType.target)

pb_sce_layer_top <- pb_sce_layer[layer_markers_top$gene, ]


row_ha <- rowAnnotation(
    cell_type = layer_markers_top$cellType.target,
    col = list(cell_type = metadata(sce)$cell_type_colors_layer)
)

column_ha <- HeatmapAnnotation(
    ncells = anno_barplot(pb_sce_layer$ncells), cell_type = pb_sce_layer$cellType_layer,
    col = list(cell_type = metadata(sce)$cell_type_colors_layer)
)

png(here(plot_dir, "heatmap_test_complex.png"), height = 800, width = 800)
Heatmap(assays(pb_sce_layer_top)$counts,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    right_annotation = row_ha,
    top_annotation = column_ha
)
dev.off()

table(sce$cellType_layer, sce$cellType_azimuth)

## plot expression of top marker for each cell type
layer_markers_best <- markers_mean_ratio |>
    filter(rank_ratio == 1) |>
    mutate(gene_anno = paste(cellType.target, "-", gene))

layer_markers_best_expression <- my_plotExpression(sce, layer_markers_best$gene, cat = "cellType_layer", fill_colors = metadata(sce)$cell_type_colors_layer)
ggsave(layer_markers_best_expression, filename = here(plot_dir, "markers_layer_best.png"))

## Add annoations
layer_df <- as.data.frame(colData(sce))[, "cellType_layer", drop = FALSE]
expression_long <- reshape2::melt(as.matrix(assays(sce)$logcounts[layer_markers_best$gene, ]))

expression_long2 <- cbind(expression_long, layer = layer_df[expression_long$Var2, ]) |>
    left_join(layer_markers_best, by = c(Var1 = "gene")) |>
    filter(!is.na(layer))

expression_long2 |> count(gene_anno)

layer_markers_best_expression2 <- expression_violin <- ggplot(data = expression_long2, aes(x = layer, y = value, fill = layer, color = layer)) +
    geom_violin(scale = "width") +
    facet_wrap(~gene_anno, nrow = 3) +
    labs(
        y = paste0("Expression logcounts")
    ) +
    theme_bw() +
    theme(
        legend.position = "None", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(face = "italic")
    ) +
    stat_summary(
        fun = median,
        geom = "crossbar",
        width = 0.3,
        color = "black"
    ) +
    scale_fill_manual(values = metadata(sce)$cell_type_colors_layer) +
    scale_color_manual(values = metadata(sce)$cell_type_colors_layer)

ggsave(layer_markers_best_expression2, filename = here(plot_dir, "markers_layer_best_anno.png"), width = 12)
