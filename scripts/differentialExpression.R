if (!requireNamespace("BiocManager", quietly = TRUE))
  installlibrary(tximport)
.packages("BiocManager")

BiocManager::install("EnrichmentBrowser")
BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("readr")

library(tximport)
library(tximportData)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(readr)
library(DESeq2)

dir <- system.file("extdata", package = "tximportData")
list.files(dir)
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
all(file.exists(files))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
head(tx2gene)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
names(txi)
head(txi$counts)
txi.tx <- tximport(files, type = "salmon", txOut = TRUE)
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)
files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample", 1:6)
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
head(txi.salmon$counts)

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
files


library(EnrichmentBrowser)
exprs.file <- "/nfs/home/students/ciora/methods/diffExp/exprs.txt"
pdat.file <- "/nfs/home/students/ciora/methods/diffExp/p_data.txt"
fdat.file <- "/nfs/home/students/ciora/methods/diffExp/f_data.txt"
de.method <- "limma"
out.file <- "/nfs/home/students/ciora/methods/diffExp/output/limma.out"

message("Reading data ...")
eset <- readSE(exprs.file, pdat.file, fdat.file)

message("DE analysis ...")
eset <- deAna(eset, de.method=de.method, padj.method="none")

de.tbl <- rowData(eset)
de.tbl <- de.tbl[, c("ENTREZID", "FC", "PVAL")]
colnames(de.tbl) <- c("GENE.ID", "log2FC", "RAW.PVAL")

write.table(de.tbl, file=out.file, row.names=FALSE, quote=FALSE, sep="\t")
message("DE table written to ", out.file)
