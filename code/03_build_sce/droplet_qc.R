library("SingleCellExperiment")
# library("DropletUtils")
library("tidyverse")
library("patchwork")
library("here")
library("sessioninfo")

## plot set up
my_theme <- theme_bw() +
  theme(text = element_text(size=15))


## Load raw data
load(here("processed-data", "sce", "sce_raw.Rdata"), verbose = TRUE)
length(table(sce$Sample))

droplet_score_fn <- list.files(here("processed-data", "03_build_sce","droplet_scores"),
                               full.names = TRUE)

names(droplet_score_fn)  <- gsub("droplet_scores_|.Rdata","",basename(droplet_score_fn))

e.out <- lapply(droplet_score_fn, function(x) get(load(x)))

## check out n empty with boxplot
FDR_cutoff <- 0.001

drop_summary <- stack(map_int(e.out, nrow)) %>%
  rename(total_n = values) %>% 
  left_join(stack(map_int(e.out, ~sum(.x$FDR < FDR_cutoff, na.rm = TRUE))) %>%
              rename(non_empty = values)) %>%
  select(Sample = ind, total_n, non_empty)

write_csv(drop_summary, file = here("processed-data", "03_build_sce","drop_summary.csv"))

drop_summary %>%
  arrange(non_empty)

summary(drop_summary$non_empty)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 2989    3682    4134    4461    5270    6269 

## compare to default drop empty
drop_summary_defualt <- read.csv(here("processed-data", "03_build_sce","drop_summary_default.csv"))

drop_summary %>% 
  left_join(drop_summary_defualt) %>% 
  mutate(d = non_empty_default - non_empty) %>%
  arrange(d)

#         Sample total_n non_empty non_empty_default     d
# 1   Br2743_ant  655424      3198              2975  -223
# 2   Br6432_ant  736110      3432              3402   -30
# 3   Br6471_ant 1148322      3622              3594   -28
# 4  Br8492_post  965917      2989              3678   689
# 5   Br2720_mid 1152093      3743              4435   692
# 6   Br6522_mid  604278      4035              4897   862
# 7   Br2743_mid 1495673      3426              4837  1411
# 8  Br6522_post  681586      4295              6066  1771
# 9   Br3942_ant 1692299      6269             52946 46677
# 10  Br6471_mid 1653348      5352             56272 50920
# 11  Br8492_mid 1302254      5189             56478 51289
# 12 Br2720_post  922885      5950             57625 51675
# 13  Br8325_ant 1855444      5951             58687 52736
# 14  Br8667_ant 1875106      5806             70117 64311
# 15  Br3942_mid 2209670      4309             74365 70056
# 16  Br8325_mid 2192531      4990             77615 72625
# 17  Br6423_ant 1521781      3999             88179 84180
# 18  Br8667_mid 2154338      4134             88660 84526
# 19 Br6423_post 1786747      4067             91276 87209

drop_barplot <- drop_summary %>%
  mutate(empty = total_n - non_empty) %>%
  select(-total_n) %>%
  pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
  ggplot(aes(x = Sample, y = n_drop, fill = drop_type))+
  geom_col() +
  scale_y_continuous(trans='log10') +
  my_theme +
  theme(axis.text.x=element_text(angle=45, hjust = 1))

ggsave(drop_barplot, filename = here("plots","03_build_sce", "drop_barplot.png"), width = 9)

## Check empty droplet results
map(e.out, ~addmargins(table(Signif = .x$FDR <= FDR_cutoff, Limited = .x$Limited, useNA = "ifany")))
map(e.out, ~addmargins(table(Signif = .x$FDR <= FDR_cutoff, useNA = "ifany")))


#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce),]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1] 36601 84756

## Compute QC metrics
sce <- addPerCellQC(
  sce,
  subsets = list(Mito = which(seqnames(sce) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

## Save for later
save(sce, file = here::here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

# sgejobs::job_single('droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript droplet_qc.R")

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

