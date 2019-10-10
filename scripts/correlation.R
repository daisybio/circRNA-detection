library(plyr)
library(dplyr)

sample_correspondence <- read.table("/data/home/students/ciora/data/mouse_brain_GSE100265/circRNA_miRNA_mRNA.tsv", header = T, stringsAsFactors = F, sep = "\t")

raw_bindSites <- read.table("/data/home/students/ciora/circRNA-detection/results/miRanda/bindsites_25%_filtered.tsv", header = T, stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]
allBindSites <- count(bindSites, Target, name="freq")

miRNA_expression_raw <- read.table("/data/home/students/ciora/circRNA-detection/results/miRNA/miRDeep2/miRNA_counts_all_samples.tsv", header = T, stringsAsFactors = F)
miRNA_expression <- miRNA_expression_raw[rowSums(miRNA_expression_raw == 0) <= 5, ]

circRNA_expression_raw <- read.table("/data/home/students/ciora/circRNA-detection/results/circExplorer2/circRNA_filtered_results.tsv",header = T, stringsAsFactors = F)
circRNA_expression <- circRNA_expression_raw[rowSums(circRNA_expression_raw == 0) <= 5, ]

correlations  <- data.frame(circRNA = character(),
                            miRNA = character(), 
                            pearson_R = numeric(), 
                            pval = numeric(),
                            stringsAsFactors=FALSE) 
for (i in 1:nrow(circRNA_expression)){
  # get coordinations of current circRNA
  chr <- as.character(circRNA_expression[i,1])
  start <- as.numeric(as.character(circRNA_expression[i,2]))
  end <- as.numeric(as.character(circRNA_expression[i,3]))
  strand <- as.character(circRNA_expression[i,4])
  circRNA <- paste(chr,":", start, "-", end, "_", strand, sep="")
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:28)]))
  colnames(circRNA_counts) <- "circRNA_counts"
  circRNA_counts$sample <- row.names(circRNA_counts)
  circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))
  
  # get miRNA binding sites of current circRNA
  #circRNA_bindSites <- count(bindSites[bindSites$Target==circRNA,], miRNA, name = "freq")
  #circRNA_bindSites <- circRNA_bindSites[order(-circRNA_bindSites$freq),]
  #miRNA <- as.character(circRNA_bindSites[2,1])
  for (j in 1:nrow(miRNA_expression)){
    miRNA <-as.character(miRNA_expression[j,1])
    
    # get sample counts for current miRNA
    miRNA_counts <- t(miRNA_expression[miRNA_expression$miRNA == miRNA,])
    miRNA_counts <- miRNA_counts[-1, ] 
    miRNA_counts <- as.data.frame(miRNA_counts)
    colnames(miRNA_counts) <- "miRNA_counts"
    miRNA_counts$sample <- row.names(miRNA_counts)
    miRNA_counts$miRNA_counts <- as.numeric(as.character(miRNA_counts$miRNA_counts))
    
    # compute circRNA expression vs. miRNA expression
    joined_counts <- merge(sample_correspondence, circRNA_counts, by.x="circRNA", by.y="sample")
    joined_counts <- merge(joined_counts, miRNA_counts, by.x="miRNA", by.y="sample")
    #ggscatter(joined_counts, x = "miRNA_counts", y = "circRNA_counts",
    #          add = "reg.line",  # Add regressin line
    #          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
    #          conf.int = TRUE) + stat_cor(method = "pearson")
    
    
    # compute circRNA-miRNA correlation for all samples
    cor_res <- cor.test(joined_counts$miRNA_counts, joined_counts$circRNA_counts,  method = "pearson", use = "complete.obs")
    R <- as.numeric(as.character(cor_res$estimate))
    pval <- as.numeric(as.character(cor_res$p.value))
    
    
    # write correlation info in correlations data frame
    correlations[nrow(correlations) + 1,] <- c(circRNA, miRNA, R, pval)
  }
  }
#write.table(correlations, file="/data/home/students/ciora/circRNA-detection/results/correlation/filtered_circRNA_miRNA_correlations.tsv", sep = "\t", quote = F, row.names = F)
correlations <- read.table("/data/home/students/ciora/circRNA-detection/results/correlation/filtered_circRNA_miRNA_correlations.tsv", stringsAsFactors = F, sep = "\t", header = T)
correlations_sign <- correlations[correlations$pval < 0.05,]

qplot(correlations$pearson_R,
      geom="histogram",
      fill=I("orange"), 
      col=I("orange"),
      alpha=I(.2),
      main="Correlation distribution between filtered circRNAs and miRNAs",
      xlab="Pearson correlation coefficient R",
      ylab="Filtered circRNA-miRNA pairs")
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/correlation/correlation_distribution.png", sep = ""), width = 8, height = 4)

qplot(correlations_sign$pearson_R,
      geom="histogram",
      fill=I("blue"), 
      col=I("blue"),
      alpha=I(.2),
      main="Significant correlation distribution (pval < 0.05)",
      xlab="Pearson correlation coefficient R",
      ylab="Filtered circRNA-miRNA pairs")

ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/correlation/correlation_distribution_sign.png", sep = ""), width = 8, height = 4)


# get pair with highest negative correlation
correlations_sign <- correlations_sign[order(correlations_sign$pearson_R),]
circRNA_min <- correlations_sign[1,1]
miRNA_min <- correlations_sign[1,2]

chr <- sapply(strsplit(as.character(circRNA_min),':'), "[", 1)
start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_min),':'), "[", 2),'-'), "[", 1))
end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_min),'-'), "[", 2),'_'), "[", 1))
strand <- sapply(strsplit(as.character(circRNA_min),'_'), "[", 2)

# get sample counts for current circRNA
circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:28)]))
colnames(circRNA_counts) <- "circRNA_counts"
circRNA_counts$sample <- row.names(circRNA_counts)
circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))

# get sample counts for current miRNA
miRNA_counts <- t(miRNA_expression[miRNA_expression$miRNA == miRNA_min,])
miRNA_counts <- miRNA_counts[-1, ] 
miRNA_counts <- as.data.frame(miRNA_counts)
colnames(miRNA_counts) <- "miRNA_counts"
miRNA_counts$sample <- row.names(miRNA_counts)
miRNA_counts$miRNA_counts <- as.numeric(as.character(miRNA_counts$miRNA_counts))

# compute circRNA expression vs. miRNA expression
joined_counts <- merge(sample_correspondence, circRNA_counts, by.x="circRNA", by.y="sample")
joined_counts <- merge(joined_counts, miRNA_counts, by.x="miRNA", by.y="sample")
ggscatter(joined_counts, x = "miRNA_counts", y = "circRNA_counts",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + stat_cor(method = "pearson") + ggtitle(paste(chr, ":", start,"-", end ," VS. ", miRNA_min, sep=""))


circRNA_bindSites <- count(bindSites[bindSites$Target==circRNA_min,], miRNA, name = "freq")
circRNA_bindSites <- circRNA_bindSites[order(-circRNA_bindSites$freq),]

# get pair with highest positive correlation
correlations_sign <- correlations_sign[order(-correlations_sign$pearson_R),]
circRNA_max <- correlations_sign[1,1]
miRNA_max <- correlations_sign[1,2]

chr <- sapply(strsplit(as.character(circRNA_max),':'), "[", 1)
start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_max),':'), "[", 2),'-'), "[", 1))
end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_max),'-'), "[", 2),'_'), "[", 1))
strand <- sapply(strsplit(as.character(circRNA_max),'_'), "[", 2)

# get sample counts for current circRNA
circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:28)]))
colnames(circRNA_counts) <- "circRNA_counts"
circRNA_counts$sample <- row.names(circRNA_counts)
circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))

# get sample counts for current miRNA
miRNA_counts <- t(miRNA_expression[miRNA_expression$miRNA == miRNA_max,])
miRNA_counts <- miRNA_counts[-1, ] 
miRNA_counts <- as.data.frame(miRNA_counts)
colnames(miRNA_counts) <- "miRNA_counts"
miRNA_counts$sample <- row.names(miRNA_counts)
miRNA_counts$miRNA_counts <- as.numeric(as.character(miRNA_counts$miRNA_counts))

# compute circRNA expression vs. miRNA expression
joined_counts <- merge(sample_correspondence, circRNA_counts, by.x="circRNA", by.y="sample")
joined_counts <- merge(joined_counts, miRNA_counts, by.x="miRNA", by.y="sample")
ggscatter(joined_counts, x = "miRNA_counts", y = "circRNA_counts",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + stat_cor(method = "pearson") + ggtitle(paste(chr, ":", start,"-", end ," VS. ", miRNA_min, sep=""))


circRNA_bindSites <- count(bindSites[bindSites$Target==circRNA_min,], miRNA, name = "freq")
circRNA_bindSites <- circRNA_bindSites[order(-circRNA_bindSites$freq),]
miRNA <- as.character(circRNA_bindSites[2,1])

