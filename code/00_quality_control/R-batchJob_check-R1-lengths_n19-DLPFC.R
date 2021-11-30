### MNT 29Nov2021 =======
  # Check Read 1 files for discrepant read lengths, as seen before in Tran-Maynard 2021 project:
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)

FASTQ.dir <- "/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/raw-data/FASTQ/"

### Read in abridged sample info (MNT generated for processing through CR)
samples.dlpfc <- read.table("/dcs04/lieber/lcolladotor/deconvolution_LIBD4030/DLPFC_snRNAseq/raw-data/sample_libs_rounds2-5.tsv",
                            sep="\t", header=F)$V1

# Add in Round1 samples (Leo had processed these & adjusted the FASTQ dir names: '-' instead of '_')
samples.dlpfc <- c(paste0(1:3, "c-k"), samples.dlpfc)

# Add in 'round0': pilot sample from cryosection:
samples.dlpfc <- c("round0", samples.dlpfc)

R1files <- data.frame(
  sampleName = unlist(sapply(samples.dlpfc, function(x){
    rep(x,length(list.files(paste0(FASTQ.dir, x), pattern="R1")))}), use.names=F),
  
  R1 = unlist(sapply(samples.dlpfc,function(x){list.files(paste0(FASTQ.dir,x),
                                                   pattern="R1")}), use.names=F)
)
dim(R1files)  # 68: [2 lanes x 4(=='round0' + 3 'Round1')] + [4 lanes x 15(rounds2-5)]

for(i in 1:nrow(R1files)){
  cat(paste0("Checking R1 length distribution for: ", R1files[i,2], "\n"))
  temp.R1s <- readFastq(paste0(FASTQ.dir, R1files[i,1], "/", R1files[i,2]),
                        withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}
sessionInfo()
