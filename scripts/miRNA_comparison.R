library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

my_folder <- "/nfs/home/students/ciora/methods/miRDeep2/output/identification/"
paper_folder <- "/nfs/home/students/ciora/data/paper_counts/miRNA/"
dataset <- read.table("/nfs/home/students/ciora/data/mouse_brain_GSE100265/miRNA_dataset", sep = "\t", header=F, stringsAsFactors = F)
my_samples <- c()
paper_samples <- c()
for(i in (1:nrow(dataset))){
  my_samples <- append(my_samples,dataset[i,1])
  paper_samples  <- append(paper_samples, dataset[i,2])
}

whole_dataset_results <- NULL
plots <- c()
for(i in 1:length(my_samples)){
  # get my results
  my_sample <- my_samples[i]
  print(i)
  current_folder <- paste(my_folder,  my_sample, '_1/', sep="")
  my_table <- read.table(dir(current_folder, full.names=T, pattern="miRNAs_expressed"), sep = "\t", header=F, stringsAsFactors = F)
  names(my_table) <- c("miRNA", "my_counts", "precursor", "total", "seq", "seq(norm)")
  my_counts <- my_table[,c(1,2)]
  # Sum up counts which come from the same miRNA because of different precursors
  my_counts <- aggregate(my_counts$my_counts, by=list(Category=my_counts$miRNA), FUN=sum)
  names(my_counts) <- c("miRNA", my_sample)
  if(is.null(whole_dataset_results)){
    whole_dataset_results <- my_counts
  } else {
    whole_dataset_results  <- merge(whole_dataset_results, my_counts, all = T, by = "miRNA")
  }
  
  # get paper's results
  paper_sample <- paper_samples[i]
  paper_counts <- read.table(dir(paper_folder, full.names=T, pattern=paper_sample), sep = "\t", header=F, stringsAsFactors = F)
  names(paper_counts) <- c("miRNA", "paper_counts")
  
  # compare my and paper's results
  names(my_counts) <- c("miRNA", "my_counts")
  comparedCounts <- merge(my_counts, paper_counts, all = T, by = "miRNA")
  comparedCounts[is.na(comparedCounts)] <- 0
  
  comparedCounts[,c(2,3)] <- comparedCounts[,c(2,3)]/1000000
   plots[[i]] <- ggscatter(comparedCounts, x = "my_counts", y = "paper_counts",
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE # Add confidence interval
   ) + stat_cor(method = "pearson") + ggtitle(paste("sample",  i)) + xlab("Current counts (millions)") + ylab("Ryback counts (millions)")
   plots[[i]]
  # ggsave(paste("/nfs/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_mapping_mine-paper-correlation_", my_sample, ".png", sep = ""), width = 8, height = 4)
  
   
   # # Consistency plots for each sample
   # ggplot(comparedCounts, aes(x=my_counts, y=paper_counts)) +
   #   geom_point() + xlim(0,2500) + ylim(0,2500) +
   #   ggtitle(paste("Consistency of miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
   #   xlab("My read counts") + ylab("Paper read counts") 
   # ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_mapping_mine_vs_paper_", my_sample, "_small.png", sep = ""), width = 8, height = 4)
  # 
  # ggplot(comparedCounts, aes(x=my_counts, y=paper_counts)) + geom_point() +
  #   ggtitle(paste("Consistency of miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
  #   xlab("My read counts") + ylab("Paper read counts")
  # ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_mapping_mine_vs_paper_", my_sample, "_whole.png", sep = ""), width = 8, height = 4)
   
   nonzero_counts <- my_counts[my_counts$my_counts > 0,]
   colnames(nonzero_counts)[2] <- "reads"
   qplot(nonzero_counts$reads,
         geom="histogram",
         fill=I("red"), 
         col=I("red"),
         alpha=I(.2),
         binwidth=100,
         main=paste("Read distribution for miRNA mapping (", my_sample, ")"),
         xlab="Read counts (> 0)",
         ylab="Number of miRNAs")
   
  #ggsave(paste("/nfs/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_read_distribution_", my_sample, "_whole.png", sep = ""), width = 8, height = 4)
  
  qplot(nonzero_counts$reads,
        geom="histogram",
        fill=I("red"), 
        col=I("red"),
        alpha=I(.2),
        binwidth=100,
        main=paste("Read distribution for miRNA mapping (", my_sample, ")"),
        xlab="Read counts (> 0)",
        ylab="Number of miRNAs",
        xlim=c(1,5000),
        ylim=c(0,100))
   #ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_read_distribution_", my_sample, "_small.png", sep = ""), width = 8, height = 4)
 
  # # read distribution mine vs paper
  # ggplot(comparedCounts, aes(x=my_counts)) + geom_histogram(fill=rgb(1,0,0,0.3)) + 
  #  xlim(0,2500) + ylim(0,90) + geom_histogram(aes(x=paper_counts), fill=rgb(0,0,1,0.3))
}
 whole_dataset_results[is.na(whole_dataset_results)]  <- 0
 #write.table(whole_dataset_results, "/data/home/students/ciora/circRNA-detection/results/miRDeep2/miRNA_counts_all_samples.tsv", quote = F, sep = "\t", row.names = F)
 
 png("/nfs/home/students/ciora/circRNA-detection/plots/miRNA/consistency_samples_1-12.png", width = 1200, height = 700)
 grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]],  layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)), top = textGrob("Consistency of miRNAs counts (paper vs. mine)", gp = gpar(fontsize = 17, font = 1)))
 dev.off()
 png("/nfs/home/students/ciora/circRNA-detection/plots/miRNA/consistency_samples_13-23.png", width = 1200, height = 700)
 grid.arrange(plots[[13]], plots[[14]], plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], plots[[21]], plots[[22]], plots[[23]], plots[[23]],  layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8), c(9, 10, 11, 12)), top = textGrob("Consistency of miRNAs counts (paper vs. mine)", gp = gpar(fontsize = 17, font = 1)))
 dev.off()
 
 
  



