## Combine our metrics with those from Tran et al

library("here")
library("sessioninfo")
library("ggpubr")

## Load previous files
load(
    here(
        "processed-data",
        "02_cellranger_metrics",
        "tran_metrics.Rdata"
    ),
    verbose = TRUE
)

load(
    here(
        "processed-data",
        "02_cellranger_metrics",
        "cellranger_metrics.Rdata"
    ),
    verbose = TRUE
)

## Combine metrics
all_metrics <- merge(tran_metrics, cellranger_metrics, all = TRUE)
dim(all_metrics)

all_metrics$Sample.ID <- basename(gsub("outs/metrics_summary.csv", "", all_metrics$metrics_csv))

## Compare metrics one at a time across studies
## Used https://rpkgs.datanovia.com/ggpubr/reference/ggboxplot.html
pdf(
    here(
        "plots",
        "02_cellranger_metrics",
        "cross_study_cellranger_metrics_boxplots.pdf"
    ),
    useDingbats = FALSE,
    width = 12,
    height = 10
)
for (i in colnames(all_metrics)[-which(colnames(all_metrics) %in% c("Sample.ID", "metrics_csv", "Q30.Bases.in.Sample.Index", "set"))]) {
    set.seed(20210831)
    p <-
        ggboxplot(
            all_metrics,
            x = "set",
            y = i,
            color = "set",
            palette = "Dark2",
            add = "jitter",
            shape = "set",
            label = "Sample.ID",
            repel = TRUE,
            font.label = list(size = 6),
            legend = "none",
            ggtheme = theme_pubr(base_size = 15)
        )
    print(p)
}
dev.off()

# sgejobs::job_single('combine_metrics', create_shell = TRUE, queue= 'bluejay', memory = '5G', command = "Rscript 03_combine_metrics.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
