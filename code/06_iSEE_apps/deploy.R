library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
options(rsconnect.max.bundle.size = 1024^3 * 4.9)
rsconnect::deployApp(
    appFiles = c("app.R", "sce_DLPFC_annotated/", "initial.R"),
    appName = "DLPFC_snRNA-seq_2022",
    account = "libd",
    server = "shinyapps.io"
)
