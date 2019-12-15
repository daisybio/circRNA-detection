library(ggplot2)

#folder <- "/nfs/home/students/ciora/methods/salmon/output/"
folder <- "/nfs/home/students/ciora/methods/salmon/sanity_check/output/"
samples <- c("SRR5144100", "SRR5144101", "SRR5144102", "SRR5144103", "SRR5144104", "SRR5144105", "SRR5144106", "SRR5144107", "SRR5144108", "SRR5144109", "SRR5144110", "SRR5144111", "SRR5144112", "SRR5144113", "SRR5144115", "SRR5144116", "SRR5144117", "SRR5144118", "SRR5144119", "SRR5144120", "SRR5144121","SRR5144122", "SRR5144123", "SRR5144124")
finaldata <- NULL
readNumber <- NULL
unnorm_counts <- NULL
for (i in 1:length(samples)){
  path <- paste(folder, samples[i], "_1/quant.sf", sep ="")
  sample_out <- read.table(path, sep = "\t", header = T)
  curr_counts <- sample_out[,c(1,5)]
  colnames(curr_counts) <- c("transcript", samples[i])
  colnames(sample_out)[c(4,5)]<-c(paste(samples[i], "_TPM", sep = ""), paste(samples[i], "_NumReads", sep = ""))
  if(is.null(finaldata)){
    finaldata <- sample_out
    readNumber <- sample_out[,4]
    unnorm_counts <- curr_counts
  } else {
    finaldata <- merge(finaldata, sample_out, by = c("Name", "Length","EffectiveLength"), all = T)
    readNumber <- append(readNumber,sample_out[,4])
    unnorm_counts <- merge(unnorm_counts, curr_counts, by = "transcript", all = T)
    }
  }
finaldata[is.na(finaldata)] <- 0
readNumber <- data.frame(readNumber)
readNumber  <- readNumber[readNumber[,1]>0,]
readNumber <- data.frame(readNumber)
#write.table(finaldata, "/nfs/home/students/ciora/methods/salmon/mRNA_results.tsv", quote = F, sep = "\t", row.names = F)
#write.table(unnorm_counts, "/nfs/home/students/ciora/methods/salmon/mRNA_unnormalized.tsv", quote = F, sep = "\t", row.names = F)
exprs <- unnorm_counts[,-1]
f_data <- data.frame(unnorm_counts[,1])
f_data$t <- unnorm_counts[,1]
#write.table(f_data, "/nfs/home/students/ciora/methods/diffExp/f_data.txt", quote = F, sep = "\t", row.names = F, col.names = F)
#write.table(exprs, "/nfs/home/students/ciora/methods/diffExp/exprs.txt", quote = F, sep = "\t", row.names = F, col.names = F)
ggplot(readNumber, aes(x=readNumber)) + geom_histogram(aes(y=..density..*100) , colour="black", fill="red") + ggtitle("Read distribution of mRNA(TPM)")+
  xlab("#reads") + ylab("density (%)") + geom_vline(aes(xintercept=5),color="blue", linetype="dashed", size=1) + xlim(0,50)
ggsave("/nfs/home/students/ciora/plots/read_distribution_mRNA.png", width = 6, height = 4)
dev.off()