#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

dataset_path <- "/nfs/home/students/ciora/testing_data/big_testing/dataset.tsv"
output_dir <- "/nfs/home/students/ciora/testing_data/big_testing/after_mapping/"

dataset <- read.table(dataset_path, sep = "\t", header=F, stringsAsFactors = F)
samples <- dataset$V1

library(plyr)
library(dplyr)
library(data.table)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)

plot_folder <- paste0(output_dir, "/results/sponging/plots/")
statistics_file <- paste0(output_dir, "/results/sponging/sponging_statistics.txt")
file.create(statistics_file)
cat("Sponging statistics",file=statistics_file,append=TRUE, sep="\n")
raw_bindSites <- read.table(paste0(output_dir, "/results/binding_sites/output/bindsites_25%_filtered.tsv"), header = T, sep = "\t", stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]
allBindSites <- count(bindSites, Target, name="freq")
bindsitDT  <- data.table(bindSites)

zero_counts_sample_cutoff <- ceiling(0.2*nrow(dataset))

miRNA_expression_raw <- read.table(paste0(output_dir, "/results/miRNA/miRNA_counts_all_samples_raw.tsv"), header = T, stringsAsFactors = F)
miRNA_expression <- miRNA_expression_raw[rowSums(miRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]

circRNA_expression_raw <- read.table(paste0(output_dir, "/results/circRNA/circRNA_counts_all_samples_filtered.tsv"),header = T, stringsAsFactors = F)
circRNA_expression <- circRNA_expression_raw[rowSums(circRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]

cat(paste0("Number of filtered miRNAs used for the correlation analysis: ", nrow(miRNA_expression)," (from ", nrow(miRNA_expression_raw), " unfiltered)"),file=statistics_file,append=TRUE, sep="\n")
cat(paste0("Number of filtered circRNAs used for the correlation analysis: ", nrow(circRNA_expression)," (from ", nrow(circRNA_expression_raw), " unfiltered)"),file=statistics_file,append=TRUE, sep="\n")

correlations <- data.frame(circRNA = character(),
                           miRNA = character(), 
                           circRNA_miRNA_ratio = numeric(),
                           miRNA_binding_sites = numeric(),
                           pearson_R = numeric(), 
                           corr_pval = numeric(),
                           RSS_norm = numeric(),
                           intercept = numeric(),
                           intercept_pval = numeric(),
                           slope = numeric(),
                           slope_pval = numeric(),
                           adj_r_squared = numeric(),
                           stringsAsFactors=FALSE) 
for (i in 1:nrow(circRNA_expression)){
  # get coordinations of current circRNA
  chr <- as.character(circRNA_expression[i,1])
  start <- as.numeric(as.character(circRNA_expression[i,2]))
  end <- as.numeric(as.character(circRNA_expression[i,3]))
  strand <- as.character(circRNA_expression[i,4])
  circRNA <- paste(chr,":", start, "-", end, "_", strand, sep="")
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:(5+length(samples)-1))]))
  colnames(circRNA_counts) <- "circRNA_counts"
  circRNA_counts$sample <- row.names(circRNA_counts)
  circRNA_counts$circRNA_counts <- as.numeric(as.character(circRNA_counts$circRNA_counts))
  
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
    joined_counts <- merge(miRNA_counts, circRNA_counts, by="sample")
    
    # analyse circRNA/miRNA ratio
    mean_circRNA_counts <- mean(joined_counts$circRNA_counts)
    mean_miRNA_counts <- mean(joined_counts$miRNA_counts)
    circRNA_miRNA_ratio <- mean_circRNA_counts/mean_miRNA_counts
    
    # compute number of miRNA binding sites on circRNA
    mirna <-miRNA
    binding_sites <- nrow(bindsitDT[miRNA == mirna & Target == circRNA])
    
    # compute circRNA-miRNA correlation for all samples
    cor_res <- cor.test(joined_counts$miRNA_counts, joined_counts$circRNA_counts,  method = "pearson", use = "complete.obs")
    corr_R <- as.numeric(as.character(cor_res$estimate))
    corr_pval <- as.numeric(as.character(cor_res$p.value))
    
    # compute linear regression
    regression_model <- lm(miRNA_counts~circRNA_counts, data = joined_counts)
    intercept <- summary(regression_model)$coefficients[1,1]
    intercept_pval <- summary(regression_model)$coefficients[1,4]
    slope <- summary(regression_model)$coefficients[2,1]
    slope_pval <- summary(regression_model)$coefficients[2,4]
    adj_r_squared <- summary(regression_model)$adj.r.squared
    
    # compute residuals sum of squares
    # normalize counts for residuals sum of squares
    normalized_counts <- joined_counts[,c("circRNA_counts", "miRNA_counts")]
    min_circRNA_counts <- min(normalized_counts$circRNA_counts)
    max_circRNA_counts <- max(normalized_counts$circRNA_counts)
    normalized_counts[,"circRNA_counts"] <- (normalized_counts[,"circRNA_counts"] - min_circRNA_counts)/(max_circRNA_counts - min_circRNA_counts)
    min_miRNA_counts <- min(normalized_counts$miRNA_counts)
    max_miRNA_counts <- max(normalized_counts$miRNA_counts)
    normalized_counts[,"miRNA_counts"] <- (normalized_counts[,"miRNA_counts"] - min_miRNA_counts)/(max_miRNA_counts - min_miRNA_counts)
    norm_reg_model <- lm(miRNA_counts~circRNA_counts, data = normalized_counts)
    RSS_norm <- sum(norm_reg_model$residuals^2)
    
    # write correlation info in correlations data frame
    correlations[nrow(correlations) + 1,] <- c(circRNA, miRNA, circRNA_miRNA_ratio, binding_sites, corr_R, corr_pval, RSS_norm, intercept, intercept_pval, slope, slope_pval, adj_r_squared)
  }
}
write.table(correlations, file=paste0(output_dir, "/results/sponging/filtered_circRNA_miRNA_correlation.tsv"), sep = "\t", quote = F, row.names = F)
correlations <- read.table(paste0(output_dir, "/results/sponging/filtered_circRNA_miRNA_correlation.tsv"), sep = "\t", stringsAsFactors = F, header = T)
correlations <- data.frame(correlations,  stringsAsFactors = F)


# plot unfiltered results

n_of_pairs <- nrow(correlations)

# create plot
p <- ggplot(correlations, aes(x=pearson_R)) +
  geom_histogram(colour=I("orange"), fill=I("orange"), alpha=I(.2)) +
  labs(title = "Correlation distribution",
       subtitle = paste(n_of_pairs,"circRNA-miRNA pairs"),
       caption ="unfiltered",
       x = "Pearson correlation coefficient R",
       y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
# add mean line
y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
p + geom_vline(xintercept = mean(correlations$pearson_R), linetype="dashed", 
               color = "black", size=0.7) + geom_text(aes(x=mean(correlations$pearson_R), label=paste(round(mean(correlations$pearson_R), digits = 3), "\n"), y = y_coord), vjust = 1.25, angle=90)

ggsave(paste0(plot_folder, "/correlation_distribution_unfiltered.png"), width = 4, height = 3)

correlations_processed <- correlations
# filter correlations (begining all combinations)
n_of_pairs_init <- nrow(correlations_processed)
print(n_of_pairs_init)
n_of_pairs <- n_of_pairs_init
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (unfiltered): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")


# filter for circRNA having mind. 1 binding site from that miRNA
bind_sites_filter = 1
correlations_processed <- correlations_processed[correlations_processed$miRNA_binding_sites >= bind_sites_filter,]
n_of_pairs <- nrow(correlations_processed)
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter,"): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")


# plot filtered results
p <- ggplot(correlations_processed, aes(x=pearson_R)) +
  geom_histogram(colour=I("purple"), fill=I("purple"), alpha=I(.2)) +
  labs(title = "Correlation distribution",
       subtitle = paste0(n_of_pairs," circRNA-miRNA pairs (", round(n_of_pairs/nrow(correlations)*100, digits=2),"%)"),
       caption = paste("Filter: binding sites >", bind_sites_filter),
       x = "Pearson correlation coefficient R",
       y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
# add mean line
y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
p + geom_vline(xintercept = mean(correlations_processed$pearson_R), linetype="dashed", 
               color = "black", size=0.7) + 
  geom_vline(xintercept = 0, color = "purple", size=0.5) +
  geom_text(aes(x=mean(correlations_processed$pearson_R), label=paste(round(mean(correlations_processed$pearson_R), digits = 3), ""), y = y_coord), vjust = 1.5, angle=90)

ggsave(paste0(plot_folder, "/corr_distr_filter_bindsites_", bind_sites_filter, ".png"), width = 4, height = 3)


# significant p-value < 0.05
pval_filter = 0.05
correlations_processed <- data.frame(correlations_processed[correlations_processed$corr_pval < pval_filter,])
n_of_pairs <- nrow(correlations_processed)
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", pval_filter , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")

# norm RSS < 1.5
maxRSS = 1.5
correlations_processed <- correlations_processed[correlations_processed$RSS_norm < maxRSS,]
n_of_pairs <- nrow(correlations_processed)
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", pval_filter , " & residual sum of squares < ", maxRSS , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")

# plot filtered results
p <- ggplot(correlations_processed, aes(x=pearson_R)) +
  geom_histogram(colour=I("green"), fill=I("green"), alpha=I(.2)) +
  labs(title = "Correlation distribution",
       subtitle = paste0(n_of_pairs," circRNA-miRNA pairs (", round(n_of_pairs/nrow(correlations)*100, digits=2),"%)"),
       caption = paste("Filter: pval <", pval_filter, ", normRSS >", maxRSS, ", bind_sites >", bind_sites_filter),
       x = "Pearson correlation coefficient R",
       y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
# add mean line
y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
p + geom_vline(xintercept = mean(correlations_processed$pearson_R), linetype="dashed", 
               color = "black", size=0.7) + 
  geom_vline(xintercept = 0, color = "green", size=0.5) +
  geom_text(aes(x=mean(correlations_processed$pearson_R), label=paste(round(mean(correlations_processed$pearson_R), digits = 3), ""), y = y_coord), vjust = 1.5, angle=90)

ggsave(paste0(plot_folder, "/corr_distr_filter_pval_",pval_filter,"_RSS_", maxRSS , "_bindsites", bind_sites_filter, ".png"), width = 4, height = 3)


# number of miRNA binding sites > 10
bind_sites_filter = 10
correlations_processed <- correlations_processed[correlations_processed$miRNA_binding_sites > bind_sites_filter,]
n_of_pairs <- nrow(correlations_processed)
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", pval_filter , " & residual sum of squares < ", maxRSS , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")



# plot filtered results
p <- ggplot(correlations_processed, aes(x=pearson_R)) +
  geom_histogram(colour=I("red"), fill=I("red"), alpha=I(.2)) +
  labs(title = "Correlation distribution",
       subtitle = paste0(n_of_pairs," circRNA-miRNA pairs (", round(n_of_pairs/nrow(correlations)*100, digits=2),"%)"),
       caption = paste("Filter: pval <", pval_filter, ", normRSS >", maxRSS, ", bind_sites >", bind_sites_filter),
       x = "Pearson correlation coefficient R",
       y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
# add mean line
y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
p + geom_vline(xintercept = mean(correlations_processed$pearson_R), linetype="dashed", 
               color = "black", size=0.7) + 
  geom_vline(xintercept = 0, color = "red", size=0.5) +
  geom_text(aes(x=mean(correlations_processed$pearson_R), label=paste(round(mean(correlations_processed$pearson_R), digits = 3), ""), y = y_coord), vjust = 1.5, angle=90)

ggsave(paste0(plot_folder, "/corr_distr_filter_pval_",pval_filter,"_RSS_", maxRSS , "_bindsites", bind_sites_filter, ".png"), width = 4, height = 3)


# number of miRNA binding sites > 50
over50 <- correlations_processed
bind_sites_filter = 50
over50 <- over50[over50$miRNA_binding_sites > bind_sites_filter,]
n_of_pairs <- nrow(over50)
print(n_of_pairs)
print(paste0(round(n_of_pairs/n_of_pairs_init*100,  digits=2), "%"))
cat(paste0("Number of circRNA-miRNA pairs (filter: number of binding sites >= ", bind_sites_filter," & correlation p-value < ", pval_filter , " & residual sum of squares < ", maxRSS , "): ", n_of_pairs, " (", round(n_of_pairs/n_of_pairs_init*100,  digits=2), " %)"), file=statistics_file, append=TRUE, sep="\n")

chr <- sapply(strsplit(as.character(over50$circRNA),':'), "[", 1)
start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(over50$circRNA),':'), "[", 2),'-'), "[", 1))
end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(over50$circRNA),'-'), "[", 2),'_'), "[", 1))
strand <- sapply(strsplit(as.character(over50$circRNA),'_'), "[", 2)
over50 <- cbind(chr, start, end, strand, over50)
over50 <- subset(over50, select=-c(circRNA, circRNA_miRNA_ratio))
write.table(over50, file=paste0(output_dir, "/results/sponging/corr_filter_pval_",pval_filter,"_RSS_", maxRSS , "_bindsites", bind_sites_filter, ".tsv"), sep = "\t", quote = F, row.names = F)

unique(over50$circRNA)

# plot filtered results
p <- ggplot(over50, aes(x=pearson_R)) +
  geom_histogram(colour=I("blue"), fill=I("blue"), alpha=I(.2)) +
  labs(title = "Correlation distribution",
       subtitle = paste0(n_of_pairs," circRNA-miRNA pairs (", round(n_of_pairs/nrow(correlations)*100, digits=2),"%)"),
       caption = paste("Filter: pval <", pval_filter, ", normRSS >", maxRSS, ", bind_sites >", bind_sites_filter),
       x = "Pearson correlation coefficient R",
       y = "circRNA-miRNA pairs") + xlim(c(-1.1,1.1))
# add mean line
y_coord <- max(ggplot_build(p)$data[[1]]$y)/2 # y coordinate for mean label
p + geom_vline(xintercept = mean(correlations_processed$pearson_R), linetype="dashed", 
               color = "black", size=0.7) + 
  geom_vline(xintercept = 0, color = "blue", size=0.5) +
  geom_text(aes(x=mean(correlations_processed$pearson_R), label=paste(round(mean(correlations_processed$pearson_R), digits = 3), ""), y = y_coord), vjust = 1.5, angle=90)

ggsave(paste0(plot_folder, "/corr_distr_filter_pval_",pval_filter,"_RSS_", maxRSS , "_bindsites", bind_sites_filter, ".png"), width = 4, height = 3)


# plot top negative correlation
correlations_sign <- correlations_processed
correlations_sign <- correlations_sign[order(correlations_sign$pearson_R),]
top_plots <- list()
for (i in 1:10){
  circRNA_min <- correlations_sign[i,1]
  miRNA_min <- correlations_sign[i,2]
  
  chr <- sapply(strsplit(as.character(circRNA_min),':'), "[", 1)
  start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_min),':'), "[", 2),'-'), "[", 1))
  end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_min),'-'), "[", 2),'_'), "[", 1))
  strand <- sapply(strsplit(as.character(circRNA_min),'_'), "[", 2)
  
  # get sample counts for current circRNA
  circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:ncol(circRNA_expression))]))
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
  joined_counts <- merge(miRNA_counts, circRNA_counts, by="sample")
  top_plots[[i]] <- ggscatter(joined_counts, y = "miRNA_counts", x = "circRNA_counts",
                              add = "reg.line",  # Add regressin line
                              add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                              conf.int = TRUE) + stat_cor(method = "pearson") +
    labs(subtitle = paste(chr, ":", start,"-", end ,"\n", miRNA_min, sep=""),
         x = "circRNA counts",
         y = "miRNA counts")
}

png(paste0(plot_folder, "/top8_neg_corr_", bind_sites_filter, "_", pval_filter, "_", maxRSS,".png"), width = 1200, height = 700)
grid.arrange(top_plots[[1]], top_plots[[2]], top_plots[[3]], top_plots[[4]], top_plots[[5]], top_plots[[6]], top_plots[[7]], top_plots[[9]],  layout_matrix = rbind(c(1, 2, 3, 4), c(5, 6, 7, 8)), top = textGrob("Top potential sponges with highest negative correlation"))
dev.off()

# get pair with highest negative correlation
correlations_sign <- correlations_processed
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
ggscatter(joined_counts, y = "miRNA_counts", x = "circRNA_counts",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + stat_cor(method = "pearson") + ggtitle(paste(chr, ":", start,"-", end ," VS. ", miRNA_min, sep=""))


circRNA_bindSites <- count(bindSites[bindSites$Target==circRNA_min,], miRNA, name = "freq")
circRNA_bindSites <- circRNA_bindSites[order(-circRNA_bindSites$freq),]

# get pair with highest positive correlation
correlations_sign <- correlations_processed
correlations_sign <- correlations_sign[order(-correlations_sign$pearson_R),]
circRNA_max <- correlations_sign[1,1]
miRNA_max <- correlations_sign[1,2]

chr <- sapply(strsplit(as.character(circRNA_max),':'), "[", 1)
start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_max),':'), "[", 2),'-'), "[", 1))
end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA_max),'-'), "[", 2),'_'), "[", 1))
strand <- sapply(strsplit(as.character(circRNA_max),'_'), "[", 2)

# get sample counts for current circRNA
circRNA_counts <- data.frame(t(circRNA_expression[circRNA_expression$chr == chr & circRNA_expression$start == start & circRNA_expression$stop == end & circRNA_expression$strand == strand,c(5:27)]))
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
joined_counts <- merge(circRNA_counts, miRNA_counts, by="sample")
ggscatter(joined_counts, y = "miRNA_counts", x = "circRNA_counts",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + stat_cor(method = "pearson") + ggtitle(paste(chr, ":", start,"-", end ," VS. ", miRNA_max, sep="")) +
  xlab("circRNA counts") + ylab("miRNA counts")

ggsave(paste0(plot_folder, "/top_pos_corr_", bind_sites_filter, "_", pval_filter, "_", maxRSS,".png"), width = 6, height = 4)

