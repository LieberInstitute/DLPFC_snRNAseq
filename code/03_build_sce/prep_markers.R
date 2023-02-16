markers.mathys.tran <- list(
    "neuron" = c("SYT1", "SNAP25", "GRIN1", "CAMK2A", "NRGN"),
    "excit_neuron" = c("SLC17A7", "SLC17A6", "SLC17A8"), # a7 coritical a6 sub-subcortical
    "inhib_neuron" = c("GAD1", "GAD2", "SLC32A1"),
    # Norepinephrine & serotonergic markers
    "neuron.NE" = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), # SLC6A3 - saw no DAT
    "neuron.5HT" = c("SLC6A4", "TPH1", "TPH2", "DDC"),
    # SERT, serotonin T (aka 5-HTT);
    "monoamine.metab" = c("COMT", "MAOA", "MAOB"),
    # MSN markers
    "MSNs.pan" = c("PPP1R1B", "BCL11B"), # "CTIP2")
    "MSNs.D1" = c("DRD1", "PDYN", "TAC1"),
    "MSNs.D2" = c("DRD2", "PENK"),
    ## Non-neuronal:
    "oligodendrocyte" = c("MBP", "MOBP", "PLP1"), # MBP can be messy
    "oligo_precursor" = c("PDGFRA", "VCAN", "CSPG4"),
    "microglia" = c("CD74", "CSF1R", "C3"),
    "astrocyte" = c("GFAP", "TNC", "AQP4", "SLC1A2"),
    "endothelial" = c("CLDN5", "FLT1", "VTN"),
    # Post-hoc from Tran-Maynard, et al. Neuron 2021
    "differn_committed_OPC" = c("SOX4", "BCAN", "GPR17", "TNS3"),
    "Tcell" = c("SKAP1", "ITK", "CD247"),
    "Mural" = c("COL1A2", "TBX18", "RBPMS"),
    "Macro" = c("CD163", "SIGLEC1", "F13A1")
)

markers.mathys.tran

save(markers.mathys.tran, file = here::here("processed-data", "03_build_sce", "markers.mathys.tran.Rdata"))

markers.Zhaung2022 <- list(
    "inhib_subtypes" = c(
        "SP8", # VIP
        "KLF5", # SST
        "LGI2", # PVALB
        "LAMP5"
    ) # LAMP5
)
