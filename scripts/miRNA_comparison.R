library(ggplot2)

my_folder <- "/data/home/students/ciora/methods/miRDeep2/output/identification/"
paper_folder <- "/data/home/students/ciora/data/paper_counts/miRNA/"
dataset <- read.table("/data/home/students/ciora/data/mouse_brain_GSE100265/miRNA_dataset", sep = "\t", header=F, stringsAsFactors = F)
my_samples <- c()
paper_samples <- c()
for(i in (1:nrow(dataset))){
  my_samples <- append(my_samples,dataset[i,1])
  paper_samples  <- append(paper_samples, dataset[i,2])
}

whole_dataset_results <- NULL
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
  
  # Consistency plots for each sample
  ggplot(comparedCounts, aes(x=my_counts, y=paper_counts)) +
    geom_point() + xlim(0,2500) + ylim(0,2500) +
    ggtitle(paste("Consistency of miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
    xlab("My read counts") + ylab("Paper read counts")
  ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_mapping_mine_vs_paper_", my_sample, "_small.png", sep = ""), width = 8, height = 4)
  #dev.off()
  
  ggplot(comparedCounts, aes(x=my_counts, y=paper_counts)) + geom_point() +
    ggtitle(paste("Consistency of miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
    xlab("My read counts") + ylab("Paper read counts")
  ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_mapping_mine_vs_paper_", my_sample, "_whole.png", sep = ""), width = 8, height = 4)
  #dev.off()
    
  ggplot(my_counts, aes(x=my_counts)) + geom_histogram(fill=rgb(1,0,0,0.3)) + 
    ggtitle(paste("Read counts distribution for miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
    xlab("Number of miRNAs") + ylab("Read counts")
  ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_read_distribution_", my_sample, "_whole.png", sep = ""), width = 8, height = 4)
  #dev.off()
  
  ggplot(my_counts, aes(x=my_counts)) + geom_histogram(fill=rgb(1,0,0,0.3)) + 
    ggtitle(paste("Read counts distribution for miRNA mapping using miRDeep2 for sample ", my_sample, "", sep = "")) + 
    xlab("Number of miRNAs") + ylab("Read counts") + xlim(0,2500) + ylim(0,100)
  ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/miRNA/samples/miRNA_read_distribution_", my_sample, "_small.png", sep = ""), width = 8, height = 4)
  #dev.off()
  
  # read distribution mine vs paper
  #ggplot(comparedCounts, aes(x=my_counts)) + geom_histogram(fill=rgb(1,0,0,0.3)) + 
  #  xlim(0,2500) + ylim(0,90) + geom_histogram(aes(x=paper_counts), fill=rgb(0,0,1,0.3))
}
 whole_dataset_results[is.na(whole_dataset_results)]  <- 0
 #write.table(whole_dataset_results, "/data/home/students/ciora/circRNA-detection/results/miRDeep2/miRNA_counts_all_samples.tsv", quote = F, sep = "\t", row.names = F)
 
 
 
 
  



