library(ggplot2)
library(grid)
library(gridExtra)

circRNAs_unfiltered <- read.table("/nfs/home/students/ciora/methods/circexplorer2_method/circRNA_results.tsv", sep="\t", header=T, stringsAsFactors = F)
dataset <- read.table("/nfs/home/students/ciora/data/mouse_brain_GSE100265/dataset_structure.txt", stringsAsFactors = F, header = T, sep = "\t")
rownames(circRNAs_unfiltered) <- paste0(circRNAs_unfiltered$chr,":",circRNAs_unfiltered$start,"-", circRNAs_unfiltered$stop,"_",circRNAs_unfiltered$strand)
circRNAs_unfiltered <- circRNAs_unfiltered[,-c(1,2,3,4)]
regions <- dataset$region
regions <- regions[-3]
colnames(circRNAs_unfiltered) <- regions


df <- t(data.frame(circRNAs_unfiltered))

pca <- prcomp(df, scale = T)

df_out <- as.data.frame(pca$x)
df_out$group <- sapply( substr(row.names(df), 1,6),"[[", 1  )

ggplot(df_out,aes(x=PC2,y=PC4,color=group ))+geom_point()+ggtitle("circRNA counts")

ggsave("/nfs/home/students/ciora/plots/analysis/PCA_circRNA_samples_2_4.png", width = 4, height = 3)









gene_counts <- read.table("/nfs/home/students/ciora/methods/salmon/genelevel_counts.tsv", sep="\t", header=T, stringsAsFactors = F)
colnames(gene_counts) <- c("CRBL_WT_1", "CRBL_WT_2", "CRBL_WT_3", "CRBL_KO_1", "CRBL_KO_2", "CRBL_KO_3", "CRTX_WT_1", "CRTX_WT_2", "CRTX_WT_3", "CRTX_KO_1", "CRTX_KO_2", "CRTX_KO_3","HIPP_WT_1", "HIPP_WT_2", "HIPP_WT_3", "HIPP_KO_1", "HIPP_KO_2", "HIPP_KO_3", "OLBU_WT_1", "OLBU_WT_2", "OLBU_WT_3", "OLBU_KO_1", "OLBU_KO_2", "OLBU_KO_3")
gene_counts <- gene_counts[rowSums(gene_counts) > 0,]
df <- t(gene_counts)
pca <- prcomp(df, scale = T)

df_out <- as.data.frame(pca$x)
df_out$group <- sapply( strsplit(as.character(row.names(df)), "_"), "[[", 1 )

ggplot(df_out,aes(x=PC3,y=PC4,color=group ))+geom_point()+ggtitle("gene counts")

ggsave("/nfs/home/students/ciora/plots/analysis/PCA_mRNA_samples_3_4.png", width = 4, height = 3)

