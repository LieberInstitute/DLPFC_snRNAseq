# Raw-data

# Sample (library prep) information

## Round 1

`Visium_LC_061621_Chromium_DLPFFC_Master.xlsx`

## Rounds 2-5

`snRNAseq_Deconvolution_Project_Master_Spreadsheet.xlsx`

```R
library(readxl)
sampleInfo <- as.data.frame(read_excel(here::here("raw-data", "sample_info", "/snRNAseq_Deconvolution_Project_Master_Spreadsheet.xlsx"), sheet=1, col_names=T))

# Pull out DLPFC samples (suffix: "_k"--i.e. Kelsey prepped these)
samples.dlpfc <- sampleInfo[grep("_k", sampleInfo[ ,"Sample #"]), 1:3]
# Add changed FASTQ file prefix...
samples.dlpfc$filePrefix <- gsub("_","-", samples.dlpfc[ ,"Sample #"])
samples.dlpfc$Tissue <- gsub(" ","_", samples.dlpfc$Tissue)
samples.dlpfc$Brain <- paste0("Br",samples.dlpfc$Brain)

write.table(samples.dlpfc, quote=F, col.names=F, row.names=F, sep="\t",
	    file=here::here("raw-data", "/sample_libs_rounds2-5.tsv"))

# Do same for LC (Matt processed: "_m")
samples.lc <- sampleInfo[grep("_m", sampleInfo[ ,"Sample #"]), 1:3]
# Add changed FASTQ file prefix...
samples.lc$filePrefix <- gsub("_","-", samples.lc[ ,"Sample #"])
samples.lc$Brain <- paste0("Br",samples.lc$Brain)

write.table(samples.lc, quote=F, col.names=F, row.names=F, sep="\t",
            file="/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/fastq/snRNA-seq/sample_libs_info.tsv")

```

# FASTQ

## Round 1


```bash
## Pilot (one sample) - calling 'round 0'
mkdir round0
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-03-15_KMay022421/Br2743_DLPFC_mid_L00*/* round0/

## 1c-k
mkdir 1c-k
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/1c_k_L00*/* 1c-k/

## 2c-k
mkdir 2c-k
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/2c_k_L00*/* 2c-k/

## 3c_k
mkdir 3c-k
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-05-21_ASpa050421/3c_k_L00*/* 3c-k/
```

## Round 2

MNT comment: See `organize_fastqs.sh` to do this iteratively down a sample 'manifest' for the remaining 15 sample libs 

and once done ->

```bash
rm organize_FASTQs.*
```

(there's definitely a more effective way to do this interactively than as an array of jobs)
