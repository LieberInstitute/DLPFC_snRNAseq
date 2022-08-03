library("SingleCellExperiment")
# library("scater")
library("tidyverse")
library("here")
library("sessioninfo")

#### Load data ####
load(here("processed-data","sce","sce_DLPFC.Rdata"), verbose = TRUE)
load(here("processed-data", "03_build_sce","cell_type_markers.Rdata"), verbose = TRUE)
# markers_1vALL
# markers_mean_ratio
# markers_mean_ratio_broad

#### set up plotting ####
source(here("code", "03_build_sce", "my_plotExpression.R"))

plot_dir <- here("plots", "05_explore_sce","03_explore_markers")
if(!dir.exists(plot_dir)) dir.create(plot_dir)

load(here("processed-data","03_build_sce","cell_type_colors.Rdata"), verbose = TRUE)
# cell_type_colors_broad
# cell_type_colors

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

(markers_1vALL_top10 <- sapply(markers_1vALL, function(x) rownames(x[1:10,])))
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

#### Plot 1vALL markers ####
markers_1vALL_top10 <- as.list(as.data.frame(markers_1vALL_top10))

my_plotMarkers(sce = sce, 
               marker_list = markers_1vALL_top10,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "markers_1vALL_top10.pdf"))

#### Mean Ratio Markers ####

markers_mean_ratio_top10 <- markers_mean_ratio |>
  group_by(cellType.target) |>
  slice(1:10) |>
  select(gene, cellType.target) |>
  unstack() |>
  as.list()


my_plotMarkers(sce = sce, 
               marker_list = markers_mean_ratio_top10,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "markers_mean_ratio_top10.pdf"))


markers_mean_ratio_broad_top10 <- markers_mean_ratio_broad |>
  group_by(cellType.target) |>
  slice(1:10) |>
  select(gene, cellType.target) |>
  unstack() |>
  as.list()

my_plotMarkers(sce = sce, 
               marker_list = markers_mean_ratio_broad_top10,
               cat = "cellType_hc",
               fill_colors = cell_type_colors,
               pdf_fn = here(plot_dir, "markers_mean_ratio_broad_top10.pdf"))

my_plotMarkers(sce = sce, 
               marker_list = markers_mean_ratio_broad_top10,
               cat = "cellType_broad_hc",
               fill_colors = cell_type_colors_broad,
               pdf_fn = here(plot_dir, "markers_mean_ratio_broad_broad_top10.pdf"))


markers_mean_ratio |> filter(gene == "MALAT1", rank_ratio <= 25)
markers_mean_ratio |> filter(gene == "MALAT1", rank_ratio <= 25)

## hm - problem cell types?
markers_mean_ratio |> filter(ratio < 1, rank_ratio <= 25) |> count(cellType.target)
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

