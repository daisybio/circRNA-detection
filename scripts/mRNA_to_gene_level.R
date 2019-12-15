library(tximportData)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

dir <- "/nfs/home/students/ciora/methods/salmon/output"
samples <- c("SRR5144100_1", "SRR5144101_1", "SRR5144102_1", "SRR5144103_1", "SRR5144104_1", "SRR5144105_1", "SRR5144106_1", "SRR5144107_1", "SRR5144108_1", "SRR5144109_1", "SRR5144110_1", "SRR5144111_1", "SRR5144112_1", "SRR5144113_1", "SRR5144115_1", "SRR5144116_1", "SRR5144117_1", "SRR5144118_1", "SRR5144119_1", "SRR5144120_1", "SRR5144121_1","SRR5144122_1", "SRR5144123_1", "SRR5144124_1")
files <- file.path(dir, samples, "quant.sf")
names(files) <- samples
all(file.exists(files))
txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
write.table(txi$counts, "/nfs/home/students/ciora/methods/salmon/genelevel_counts.tsv", quote = F, sep = "\t", row.names = T)
