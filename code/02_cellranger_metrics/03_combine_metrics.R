## Combine our metrics with those from Tran et al

library("here")
library("sessioninfo")

## Load previous files
load(here(
    "processed-data",
    "02_cellranger_metrics",
    "tran_metrics.Rdata"
),
    verbose = TRUE)

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

## Compare metrics one at a time across studies
## Used https://rpkgs.datanovia.com/ggpubr/reference/ggboxplot.html
pdf(
    here(
        "processed-data",
        "02_cellranger_metrics",
        "cross_study_cellranger_metrics_boxplots.pdf"
    ),
    useDingbats = FALSE,
    width = 10
)
for (i in colnames(all_metrics)[-c(1, ncol(all_metrics))]) {
    set.seed(20210714)
    p <-
        ggboxplot(
            all_metrics,
            x = "study",
            y = i,
            color = "study",
            palette = "Dark2",
            add = "jitter",
            shape = "study",
            label = "Sample.ID",
            repel = TRUE,
            font.label = list(size = 5),
            legend = "none",
            ggtheme = theme_pubr(base_size = 30)
        )
    print(p)
}
dev.off()



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
