library(plyr)
library(dplyr)
file <- "/data/home/students/ciora/methods/miRanda/output/circRNA_bind_sites_results.txt"
plot_folder <- "/data/home/students/ciora/circRNA-detection/plots/binding_sites/"
bindSites <- read.table(file, header = T, stringsAsFactors = F)
bindSites <- bindSites[,c(1,2)]

allBindSites <- count(bindSites, Target, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- count(distinct, Target, name="freq")

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
  