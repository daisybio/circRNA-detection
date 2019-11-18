#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied", call.=FALSE)
}
output_dir = args[1]

library(ggplot2)
library(plyr)
library(dplyr)

#output_dir <- "/nfs/home/students/ciora/testing_data/small_testing/output_whole/"

file <- paste0(output_dir, "/results/binding_sites/output/circRNA_bind_sites_results.txt")
plot_folder <- paste0(output_dir, "results/binding_sites/plots/")
raw_bindSites <- read.table(file, header = T, sep = '\t', stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]

allBindSites <- count(bindSites, Target, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- count(distinct, Target, name="freq")

write.table(raw_bindSites, file = paste0(output_dir, "/results/binding_sites/output/bindsites_raw.tsv"), sep = "\t", quote = F, row.names = F)
write.table(allBindSites, file = paste0(output_dir, "/results/binding_sites/output/bindsites_per_circRNA.tsv"), sep = "\t", quote = F, row.names = F)

ggplot(distinctBindSites, aes(x=Target, y=freq)) + geom_bar(stat="identity", color = "red", fill= "red") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle("Number of distinct binding sites per circRNA") +
  xlab(paste("", nrow(distinctBindSites), " target circRNA", sep="")) +
  ylab("Number of binding sites")
ggsave(paste(plot_folder, "distinct_bind_sites.png", sep=""), width = 6, height = 4)
dev.off()

ggplot(allBindSites, aes(x=Target, y=freq)) + geom_bar(stat="identity", color = "red", fill= "red") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  ggtitle("Number of binding sites per circRNA") +
  xlab(paste(nrow(allBindSites), " target circRNA", sep="")) +
  ylab("Number of binding sites")
ggsave(paste(plot_folder, "all_bind_sites.png", sep=""), width = 6, height = 4)
dev.off()

ggplot(allBindSites) + geom_histogram(colour="blue", fill="lightblue", aes(x=freq)) +
  ggtitle("Distribution of miRNA binding sites on circRNA targets") +
  xlab("Number of binding sites") + ylab("Number of circRNAs")
ggsave(paste(plot_folder, "all_bind_sites_histogram.png", sep=""), width = 6, height = 4)
dev.off()

ggplot(distinctBindSites) + geom_histogram(colour="green", fill="lightgreen", aes(x=freq)) +
  ggtitle("Distribution of distinct miRNA binding sites on circRNA targets") +
  xlab("Number of distinct binding sites") + ylab("Number of circRNAs")
ggsave(paste(plot_folder, "distinct_bind_sites_histogram.png", sep=""), width = 6, height = 4)
dev.off()

# miRanda score distribution
scores <- raw_bindSites[,c(1,3)]

# fiter 25% worst scores
filtered_scores <- raw_bindSites[raw_bindSites$Score > quantile(raw_bindSites$Score, 0.25),]
write.table(filtered_scores, file = paste0(output_dir, "/results/binding_sites/output/bindsites_25%_filtered.tsv"), sep = "\t", quote = F, row.names = F)
filtered_scores <- filtered_scores[,c(1,2)]

# redo plots after filtering
allBindSites <- count(filtered_scores, Target, name="freq")
distinct <- distinct(filtered_scores)
distinctBindSites <- count(distinct, Target, name="freq")

ggplot(distinctBindSites, aes(x=Target, y=freq)) + geom_bar(stat="identity", color = "red", fill= "red") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle("Number of distinct binding sites per circRNA (filtered)") +
  xlab(paste("", nrow(distinctBindSites), " target circRNA", sep="")) +
  ylab("Number of binding sites")
ggsave(paste(plot_folder, "distinct_bind_sites_filtered.png", sep=""), width = 6, height = 4)

ggplot(allBindSites, aes(x=Target, y=freq)) + geom_bar(stat="identity", color = "red", fill= "red") + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  ggtitle("Number of binding sites per circRNA (filtered)") +
  xlab(paste(nrow(allBindSites), " target circRNA", sep="")) +
  ylab("Number of binding sites")
ggsave(paste(plot_folder, "all_bind_sites_filtered.png", sep=""), width = 6, height = 4)

ggplot(allBindSites) + geom_histogram(colour="blue", fill="lightblue", aes(x=freq)) +
  ggtitle("Distribution of miRNA binding sites on circRNA targets (filtered)") +
  xlab("Number of binding sites") + ylab("Number of circRNAs")
ggsave(paste(plot_folder, "all_bind_sites_histogram_filtered.png", sep=""), width = 6, height = 4)

ggplot(distinctBindSites) + geom_histogram(colour="green", fill="lightgreen", aes(x=freq)) +
  ggtitle("Distribution of distinct miRNA binding sites on circRNA targets (filtered)") +
  xlab("Number of distinct binding sites") + ylab("Number of circRNAs")
ggsave(paste(plot_folder, "distinct_bind_sites_histogram_filtered.png", sep=""), width = 6, height = 4)


