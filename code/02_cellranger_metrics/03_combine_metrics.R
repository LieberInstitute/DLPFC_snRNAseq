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

## Subset to current samples and explore
current <- subset(all_metrics, set %in% paste0("round", 0:5))
dim(current)
# [1] 19 23
summary(current)
 # Estimated.Number.of.Cells Mean.Reads.per.Cell Median.Genes.per.Cell Number.of.Reads     Valid.Barcodes
 # Min.   :3156              Min.   : 38062      Min.   : 985          Min.   :162821623   Min.   :92.00
 # 1st Qu.:3906              1st Qu.: 47784      1st Qu.:1434          1st Qu.:275981256   1st Qu.:93.80
 # Median :5870              Median : 55811      Median :1959          Median :323114123   Median :94.70
 # Mean   :5583              Mean   : 65253      Mean   :2415          Mean   :351032590   Mean   :95.02
 # 3rd Qu.:6902              3rd Qu.: 75312      3rd Qu.:3259          3rd Qu.:483252608   3rd Qu.:96.50
 # Max.   :7848              Max.   :136799      Max.   :5146          Max.   :541289607   Max.   :97.40
 #
 # Sequencing.Saturation Q30.Bases.in.Barcode Q30.Bases.in.RNA.Read Q30.Bases.in.UMI Reads.Mapped.to.Genome
 # Min.   :21.70         Min.   :93.20        Min.   :91.80         Min.   :93.60    Min.   :91.00
 # 1st Qu.:41.35         1st Qu.:96.90        1st Qu.:94.75         1st Qu.:96.80    1st Qu.:96.70
 # Median :54.60         Median :97.00        Median :95.00         Median :96.90    Median :97.00
 # Mean   :57.78         Mean   :96.78        Mean   :94.84         Mean   :96.75    Mean   :96.62
 # 3rd Qu.:76.25         3rd Qu.:97.10        3rd Qu.:95.30         3rd Qu.:97.00    3rd Qu.:97.20
 # Max.   :88.20         Max.   :97.20        Max.   :95.60         Max.   :97.10    Max.   :97.30
 #
 # Reads.Mapped.Confidently.to.Genome Reads.Mapped.Confidently.to.Intergenic.Regions
 # Min.   :67.40                      Min.   :4.700
 # 1st Qu.:80.70                      1st Qu.:5.900
 # Median :83.40                      Median :6.100
 # Mean   :83.35                      Mean   :6.358
 # 3rd Qu.:87.35                      3rd Qu.:6.700
 # Max.   :92.20                      Max.   :8.200
 #
 # Reads.Mapped.Confidently.to.Intronic.Regions Reads.Mapped.Confidently.to.Exonic.Regions
 # Min.   :21.40                                Min.   :22.50
 # 1st Qu.:38.85                                1st Qu.:28.05
 # Median :46.00                                Median :33.00
 # Mean   :43.66                                Mean   :33.32
 # 3rd Qu.:52.40                                3rd Qu.:37.15
 # Max.   :61.00                                Max.   :49.30
 #
 # Reads.Mapped.Confidently.to.Transcriptome Reads.Mapped.Antisense.to.Gene Fraction.Reads.in.Cells Total.Genes.Detected
 # Min.   :48.20                             Min.   : 4.90                  Min.   :37.80           Min.   :28060
 # 1st Qu.:57.25                             1st Qu.: 9.55                  1st Qu.:53.70           1st Qu.:30303
 # Median :63.50                             Median :14.70                  Median :70.30           Median :31734
 # Mean   :62.03                             Mean   :14.19                  Mean   :66.56           Mean   :31241
 # 3rd Qu.:64.60                             3rd Qu.:17.60                  3rd Qu.:79.40           3rd Qu.:32496
 # Max.   :76.00                             Max.   :28.00                  Max.   :86.20           Max.   :32967
 #
 # Median.UMI.Counts.per.Cell metrics_csv            set            Q30.Bases.in.Sample.Index  Sample.ID
 # Min.   : 1323              Length:19          Length:19          Min.   : NA               Length:19
 # 1st Qu.: 2144              Class :character   Class :character   1st Qu.: NA               Class :character
 # Median : 3334              Mode  :character   Mode  :character   Median : NA               Mode  :character
 # Mean   : 5253                                                    Mean   :NaN
 # 3rd Qu.: 7366                                                    3rd Qu.: NA
 # Max.   :14627                                                    Max.   : NA
 #                                                                  NA's   :19

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
