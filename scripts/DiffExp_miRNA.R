library(DESeq2)
library(ggplot2)
library(gridExtra)
library(grid)

log2FClevel = 1
padjlevel=0.05

# preprocessing counts table
counts <- read.table("/nfs/home/students/ciora/circRNA-detection/results/miRNA/miRDeep2/miRNA_counts_all_samples.tsv", sep = "\t", header=T, stringsAsFactors = F)
rownames(counts) <- counts$miRNA
counts$miRNA <- NULL

# preprocessing dataset information
dataset_structure <- read.table("/nfs/home/students/ciora/data/mouse_brain_GSE100265/miRNA_dataset", sep = "\t", header=F, stringsAsFactors = F)
dataset <- dataset_structure[,c(1,3)]
colnames(dataset) <- c("sample","region")
dataset$condition <- sapply(strsplit(as.character(dataset_structure$V4),'_'), "[", 1)
dataset$rep <- sapply(strsplit(as.character(dataset_structure$V4),'_'), "[", 2)
coldata <- dataset[,c("condition", "region")]
rownames(coldata) <- dataset$sample
all(rownames(coldata) == colnames(counts))

#computing DESeq for WT vs. KO for every region and write results tables
regions <- unique(dataset$region)
significant_miRNAs <- data.frame(matrix(ncol = length(regions), nrow = 0))
colnames(significant_miRNAs) <- regions
for( i in (1:length(regions))){
  # DESeq for specific region
  specific_region <- regions[i]
  specific_coldata <- coldata[coldata$region == specific_region,]
  specific_counts <- counts[,rownames(specific_coldata)]
  dds <- DESeqDataSetFromMatrix(countData = specific_counts, colData = specific_coldata, design = ~ condition)
  dds <- DESeq(dds)
  resultsNames(dds)
  res <- results(dds, name="condition_WT_vs_KO")
  res
  res <- results(dds, contrast=c("condition","WT","KO"))
  res
  resOrdered <- res[order(res$pvalue),]
  summary(res)
  resOrdered
  #plotMA(res, ylim=c(-5,5))
  #idx <- identify(res$baseMean, res$log2FoldChange)
  #rownames(res)[idx]
  
  #write all results
  write.table(as.data.frame(resOrdered), quote = F, sep = "\t",
              file=paste("/nfs/home/students/ciora/circRNA-detection/results/miRNA/DiffExp/DE_miRNA_WTvsKO_", specific_region, "_all.tsv", sep = ""))
  
  #write only significant results
  resSig <- subset(resOrdered, padj < padjlevel & abs(log2FoldChange) >= log2FClevel)
  resSig <- as.data.frame(resSig)
  write.table(as.data.frame(resSig), quote = F, sep = "\t",
              file=paste("/nfs/home/students/ciora/circRNA-detection/results/miRNA/DiffExp/DE_miRNA_WTvsKO_", specific_region, "_significant_padj<", padjlevel, "log2FC>=" ,log2FClevel,".tsv", sep = ""))

  #summarize all significant miRNAs to one table
  if(nrow(resSig) > 0){
    for(k in (1:nrow(resSig))){
      miRNA <- row.names(resSig)[k]
      significant_miRNAs[miRNA,specific_region] <- miRNA
    }
  }
}

write.table(as.data.frame(significant_miRNAs), quote = F, row.names = F, sep = "\t",
            file=paste("/nfs/home/students/ciora/circRNA-detection/results/miRNA/DiffExp/significant_miRNAs_all_regions_padj<", padjlevel, "log2FC>=" ,log2FClevel, ".tsv", sep = ""))

#plotMA(res, ylim=c(-5,5))

# volcano plots
plots <- list()
for(j in 1:length(regions)){
  reg <- regions[j]
  path <- paste("/nfs/home/students/ciora/circRNA-detection/results/miRNA/DiffExp/DE_miRNA_WTvsKO_", reg, "_all.tsv", sep = "")
  table <- read.table(path, header = T)
  table$Legend <- "non-significant"
  table$Legend[abs(table$log2FoldChange) > log2FClevel & table$padj <= padjlevel] <- "significant"
  title <- paste("", reg, " region", sep = "")
  plots [[j]] <- ggplot(table, aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(aes(col=Legend), cex = 1) + 
    scale_color_manual(values = c("black", "red")) + 
    ggtitle(title)+ theme(legend.position="none") + 
    geom_hline(yintercept=-log10(padjlevel), linetype="dashed", 
                    color = "red", size=0.5) + 
    geom_vline(xintercept=log2FClevel, linetype="dashed", 
               color = "red", size=0.5) +
    geom_vline(xintercept=-log2FClevel, linetype="dashed", 
               color = "red", size=0.5)
    #+ annotate("text", label = paste("padj <=", padjlevel), 
    #         x = 0, y = -log10(padjlevel), color = "black", size = 3) +
    #annotate("text", label = paste("log2FC >", log2FClevel), 
    #         x = 1, y = 2, color = "black", size = 3, angle = 90)
  }
png(paste("/nfs/home/students/ciora/plots/miRNA/DE_miRNAs_regions_padj<", padjlevel, "log2FC>=" ,log2FClevel,".png", sep = ""), width = 700, height = 600)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], layout_matrix = rbind(c(1, 2),c(3, 4)), top = textGrob("Differentially expressed miRNAs (WT vs. KO)", gp = gpar(fontsize = 17, font = 1)))
dev.off()

