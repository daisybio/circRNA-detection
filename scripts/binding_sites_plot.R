library(plyr)
library(dplyr)
file <- "/data/home/students/ciora/methods/miRanda/output/circRNA_bind_sites_results.txt"
plot_folder <- "/data/home/students/ciora/circRNA-detection/plots/binding_sites/"
raw_bindSites <- read.table(file, header = T, stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]

allBindSites <- count(bindSites, Target, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- count(distinct, Target, name="freq")

write.table(raw_bindSites, file = "/data/home/students/ciora/circRNA-detection/results/miRanda/bindsites_raw.tsv", sep = "\t", quote = F, row.names = F)
write.table(allBindSites, file = "/data/home/students/ciora/circRNA-detection/results/miRanda/bindsites_per_circRNA.tsv", sep = "\t", quote = F, row.names = F)

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
ggplot(data=scores, aes(x=Score)) + geom_histogram(fill=I("red"), col=I("red"), alpha=I(.2)) +
  ggtitle("miRanda score distribution") + xlab("Score") + ylab("Number of binding sites") + 
  geom_vline(xintercept=quantile(scores$Score, 0.25), linetype="dashed", 
             color = "blue", size=0.5) +
  annotate("text", x = quantile(scores$Score, 0.25) - 2, y = 1000000, label = "25%", color = "blue", angle = 90)
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/binding_sites/miRanda_score_distribution.png", sep = ""), width = 8, height = 4)

# fiter 25% worst scores
filtered_scores <- raw_bindSites[raw_bindSites$Score > quantile(raw_bindSites$Score, 0.25),]
write.table(filtered_scores, file = "/data/home/students/ciora/circRNA-detection/results/miRanda/bindsites_25%_filtered.tsv", sep = "\t", quote = F, row.names = F)
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
  

  
  # compute length of every circRNA
  f = function(x, output) {
    # x is the row of type Character
    circRNA = as.character(x[1])
    start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA),':'), "[", 2),'-'), "[", 1))
    end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA),'-'), "[", 2),'_'), "[", 1))
    length <- end - start
    freq = x[2]
    print(length)
    #your code to process x
  }
  
  result <- apply(allBindSites, 1, f)
  allBindSites$length <- result
  
  ggscatter(allBindSites, x = "freq", y = "length",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE   ) + stat_cor(method = "pearson") + 
    ggtitle("Correlation: length vs. binding sites of circRNA") +
    xlab("Number of binding sites") + ylab("Length of circRNA")
  
  ggsave("/data/home/students/ciora/circRNA-detection/plots/binding_sites/corr_length_bindsites_filtered.png", width = 6, height = 4)
  
  
  
  