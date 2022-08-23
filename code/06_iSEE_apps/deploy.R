library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_DLPFC_annotated/", "initial.R"),
    appName = "DLPFC_snRNA-seq_2022",
    account = "libd",
    server = "shinyapps.io"
)

