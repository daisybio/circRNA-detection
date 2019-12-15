library(plyr)
library(dplyr)
library(Gviz)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(IRanges)


# read ATXN1 known circRNAs from circBase
circBase_ATXN1 <- read.table("circRNAs_ATXN1_mm10.tsv", stringsAsFactors = F)
circBase_ATXN1$chromosome <- sapply(strsplit(as.character(circBase_ATXN1[,1]),':'), "[", 1)
circBase_ATXN1$start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circBase_ATXN1[,1]),':'), "[", 2),'-'), "[", 1))
circBase_ATXN1$end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circBase_ATXN1[,1]),'-'), "[", 2),'_'), "[", 1))
circBase_ATXN1$strand <- circBase_ATXN1$V2
circBase_ATXN1 <- circBase_ATXN1[,-c(1,2)]

min_start <- min(circBase_ATXN1$start)
max_end <- max(circBase_ATXN1$end)

# read circRNAs detected by CircExplorer2 on the same region
circRNAs <- read.table("circRNA_filtered_results.tsv", sep = "\t", stringsAsFactors = F, header = T)
circRNA_detected <- circRNAs[,c(1:4)]
colnames(circRNA_detected) <- c("chromosome", "start", "end", "strand")
circRNA_detected <- circRNA_detected[circRNA_detected$chromosome =="chr13" & circRNA_detected$start >= min_start & circRNA_detected$end <= max_end,]

# plot known vs. detected circRNAs on ATXN1 region
chr <- "chr13"
gen <- "mm10"
atrack_circ_known <- AnnotationTrack(circBase_ATXN1, name="known circRNAs (circBase)", fill = "green3")
atrack_circ_detected <- AnnotationTrack(circRNA_detected, name="detected circRNAs", fill = "orange")
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome=gen, chromosome=chr)
plotTracks(list(itrack, gtrack, atrack_circ_known, atrack_circ_detected))

# plot expression of atxn1 circRNAs in regions and samples

circRNA_expression_raw <- read.table("circRNA_counts_all_samples_filtered.tsv",header = T, stringsAsFactors = F)
circRNA_expression <- circRNA_expression_raw[rowSums(circRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]
atxn1_expression <- circRNA_expression_raw[circRNA_expression_raw$chr == chr,]
atxn1_expression <- atxn1_expression[atxn1_expression$start >= min_start & atxn1_expression$stop <= max_end,]
rownames(atxn1_expression) <- paste0(atxn1_expression$chr, ":", atxn1_expression$start, "-", atxn1_expression$stop, "_", atxn1_expression$strand)
atxn1_expression <- atxn1_expression[,-c(1,2,3,4)]
samples <- colnames(atxn1_expression)
atxn1_expression$circRNA <- rownames(atxn1_expression)

df <- melt(atxn1_expression, id.vars = c("circRNA"),
     measure.vars = samples)
colnames(df) <- c("circRNA", "sample", "counts")
df$region <- sapply(strsplit(as.character(df$sample),'_'), "[", 1)
df$type <- sapply(strsplit(as.character(df$sample),'_'), "[", 2)
df$rep <- sapply(strsplit(as.character(df$sample),'_'), "[", 3)

WT <- df[df$type == "WT",c("circRNA", "counts", "region", "rep")]
KO <- df[df$type == "KO",c("circRNA", "counts", "region", "rep")]

ggplot(WT) + aes(x=circRNA, y=counts, fill = circRNA) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_bar(stat="identity", width = 0.7) +  
  ggtitle("Expression of detected circRNAs on ATXN1 gene (WT)") +
  xlab("sample") + ylab("counts") + facet_grid(rep~region)

ggplot(KO) + aes(x=circRNA, y=counts, fill = circRNA) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_bar(stat="identity", width = 0.7) +  
  ggtitle("Expression of circRNAs detected on ATXN1 gene (KO)") +
  xlab("sample") + ylab("counts") + facet_grid(rep~region)
