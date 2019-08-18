library(ggplot2)

folder <- "/data/home/students/ciora/methods/circexplorer2_method/final_test/"
samples <- c("SRR5720118", "SRR5720119", "SRR5720120", "SRR5720121", "SRR5720122", "SRR5720123", "SRR5720124", "SRR5720125", "SRR5720126", "SRR5720127", "SRR5720128", "SRR5720129", "SRR5720130", "SRR5720131", "SRR5720132", "SRR5720133", "SRR5720134", "SRR5720135", "SRR5720136", "SRR5720137", "SRR5720138", "SRR5720139", "SRR5720140", "SRR5720141")
alldata <- NULL
finaldata <- NULL
counts <- c()
counts_cutoff <- c()
cutoff <- 5
for (i in 1:length(samples)){
  path <- paste(folder, samples[i], "_1/circ_out/circularRNA_known.txt", sep ="")
  sample_out <- read.table(path, sep = "\t")
  counts[i] <- nrow(sample_out)
  names(counts)[i] <- samples[i]
  counts_cutoff[i] <- nrow(sample_out[sample_out[,13]>cutoff,])
  names(counts_cutoff)[i] <- samples[i]
  expression <- sample_out[,c(1,2,3,6,13,15,16)]
  colnames(expression) <- c("chr", "start", "stop", "strand", samples[i], paste("gene_", samples[i], sep = ""), paste("isoform_", samples[i], sep = ""))
  expression <- expression[,c(1,2,3,4,5)]
  if(is.null(finaldata)){
    finaldata <- expression
    alldata <- expression[,samples[i]]
  } else {
    finaldata <- merge(finaldata, expression, by = c("chr", "start", "stop", "strand"), all = T)
    alldata <- append(alldata, expression[,samples[i]])
    }
  
}
finaldata[is.na(finaldata)] <- 0
#write.table(finaldata, "/data/home/students/ciora/methods/circexplorer2_method/circRNA_results.tsv", quote = F, sep = "\t", row.names = F)

counts <- data.frame(sample=names(counts),counts)
counts_cutoff <- data.frame(sample=names(counts_cutoff), counts_cutoff)
dataset <- read.table("/data/home/students/ciora/data/mouse_brain_GSE100265/dataset_structure.txt", sep = "\t", header=T)
sub_dataset <- merge(counts, dataset, by.x = "sample", by.y="SRA", all= T)
sub_dataset[is.na(sub_dataset)] <- 0
sub_dataset_cutoff <- merge(counts_cutoff, dataset, by.x = "sample", by.y="SRA", all= T)
sub_dataset_cutoff[is.na(sub_dataset_cutoff)] <- 0

#prep_dataset <- data.frame(matrix(0, nrow = 6, ncol = 4))
#colnames(prep_dataset) <- c("cerebellum", "cortex", "hippocampus", "olfactory bulb")
#rownames(prep_dataset) <- c("KO_rep1", "KO_rep2", "KO_rep3", "WT_rep1", "WT_rep2", "WT_rep3")
#for(i in 1:nrow(sub_dataset)){
#  region <- as.character(sub_dataset[i,"region"])
#  type <- paste(as.character(sub_dataset[i,"type"]), as.character(sub_dataset[i,"rep"]), sep="_")
#  prep_dataset[type,region] <- sub_dataset[i,"counts"]
#}

# diversity before cutoff
ggplot(sub_dataset) + aes(x=type_rep, y=counts) + 
  geom_bar(stat="identity", width = 0.7, fill = "salmon") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Diversity of circular RNAs") +
  xlab("Sample") + ylab("# circular RNA") + facet_grid(~region)
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/amountCircRNA_before_", cutoff, ".png", sep = ""), width = 8, height = 4)
dev.off()

#diversity after cutoff
ggplot(sub_dataset_cutoff) + aes(x=type_rep, y=counts_cutoff) + 
  geom_bar(stat="identity", width = 0.7, fill="salmon") + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle(paste("Diversity of circular RNAs with > ", cutoff, " reads", sep = "")) +
  xlab("Sample") + ylab("# circular RNA") + facet_grid(~region)
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/amountCircRNA_after_", cutoff, ".png", sep = ""), width = 8, height = 4)
dev.off()

# diversity before and after
combined_dataset <- merge(sub_dataset[ , c("sample", "counts", "region", "type_rep")], sub_dataset_cutoff[ ,c("sample", "counts_cutoff", "region", "type_rep")], by = c("sample", "region","type_rep"), )
combined_dataset$difference <- combined_dataset[,"counts"] - combined_dataset[,"counts_cutoff"]
stacked_data <- NULL
stacked_data <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(stacked_data) <- c("sample",  "counts", "region", "type_rep", "cutoff")
for (i in 1:nrow(combined_dataset)){
  row <- nrow(stacked_data) + 1
  stacked_data[c(row,row+1), "sample"] <- as.character(combined_dataset[i, "sample"])
  stacked_data[c(row, row+1), "region"] <- as.character(combined_dataset[i, "region"])
  stacked_data[c(row, row+1), "type_rep"] <- as.character(combined_dataset[i, "type_rep"])
  stacked_data[row, "counts"] <- combined_dataset[i, "difference"]
  stacked_data[row+1, "counts"] <- combined_dataset[i, "counts_cutoff"]
  stacked_data[row, "cutoff"] <- paste ("<= ", cutoff, " reads", sep = "" )
  stacked_data[row+1, "cutoff"] <- paste ("> ", cutoff, " reads", sep = "" )
}
ggplot(stacked_data) + aes(x=type_rep, y=counts, fill = cutoff)+
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=c("#A0A0A0", "salmon")) +
  ggtitle("Diversity of circular RNAs before and after cutoff") +
  xlab("Sample") + ylab("# circular RNA") + facet_grid(~region)
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/amountCircRNA_before_after_", cutoff, ".png", sep = ""), width = 8, height = 4)
dev.off()

# read distribution before cutoff
alldata  <- data.frame(alldata)
ggplot(alldata, aes(x=alldata)) + geom_histogram(aes(y=..density..*100) , colour="black", fill="yellow", binwidth = 1, boundary = 0.5) +
  xlim(0,50) + ggtitle("Read distribution before cutoff (circRNA, all samples)") +
  xlab("Number of reads") + ylab("density (%)") + 
  geom_vline(aes(xintercept=cutoff), color="blue", linetype="dashed", size=1) +
  annotate("text", x = cutoff + 1, y = 30, label = paste("cutoff = ", cutoff, sep = ""), 
           color = "blue", angle = 270)
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/read_distribution_before_cutoff_", cutoff, ".png", sep=""), width = 6, height = 4)
dev.off()

# read distribution after cutoff
cutoffdata <- data.frame(alldata[alldata[,1] > cutoff ,])
colnames(cutoffdata)[1] <- "reads"
ggplot(cutoffdata, aes(x=reads)) + geom_histogram(aes(y=..density..*100),colour="black", fill="yellow", binwidth = 1)+
  xlim(cutoff,50) + ggtitle(paste("Read distribution after cutoff > ", cutoff, " (circRNA, all samples)", sep = "")) +
  xlab("Number of reads") + ylab("density (%)")
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/read_distribution_after_cutoff_", cutoff, ".png", sep=""), width = 6, height = 4)
dev.off()

t <- data.frame(table(cutoffdata))

