library(VennDiagram)

mm9_known_circRNA <- read.table("/data/home/students/ciora/circRNA-detection/data/mm9_known_circRNA.txt", sep = "\t", header = F, stringsAsFactors = F)
#liftOver_input <- paste(mm9_circRNA$V1, ":", mm9_circRNA$V2, "-", mm9_circRNA$V3, sep = "")
#write.table(liftOver_input ,"/data/home/students/ciora/circRNA-detection/data/mm9_circRNA_processed.txt", quote = F, row.names = F, col.names = F)
mm10_known_circRNA <- read.table("/data/home/students/ciora/circRNA-detection/data/mm10_known_circRNA.bed", sep = "\t", header = F, stringsAsFactors = F)

# reference circBase all circRNAs known (doesn't contain Rybak)
mm10_known_circRNA$chr <- sapply(strsplit(as.character(mm10_known_circRNA$V1),':'), "[", 1)
mm10_known_circRNA$start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(mm10_known_circRNA$V1),':'), "[", 2),'-'), "[", 1))
mm10_known_circRNA$end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(mm10_known_circRNA$V1),':'), "[", 2),'-'), "[", 2))
mm10_known_circRNA$V1 <- NULL

#reference circBase Rybak experiment
mm10_Rybak_circRNA <- read.table("/data/home/students/ciora/circRNA-detection/data/mm10_Rybak.bed", sep = "\t", header = F, stringsAsFactors = F)[,c(1,2,3)]
colnames(mm10_Rybak_circRNA) <- c("chr", "start", "end")

mm10_merged_known_circRNA <- merge(mm10_Rybak_circRNA, mm10_known_circRNA)

# circRNAs detected by me using circExplorer2
mm10_my_circRNA <- read.table("/data/home/students/ciora/circRNA-detection/results/circExplorer2/circRNA_filtered_results.tsv", header = T, stringsAsFactors = F, sep = "\t")

result <- mm10_my_circRNA[,c(1,2,3)]
for(i in (1:nrow(mm10_my_circRNA))){
  my_chr <- mm10_my_circRNA[i,1]
  my_start <- as.numeric(mm10_my_circRNA[i,2])
  my_end <- as.numeric(mm10_my_circRNA[i,3])
  sub_Rybak <- mm10_Rybak_circRNA[(mm10_Rybak_circRNA[,"chr"]==my_chr & mm10_Rybak_circRNA[,"start"]>=my_start & mm10_Rybak_circRNA[,"end"]<=my_end),]
  sub_other <- mm10_known_circRNA[(mm10_known_circRNA[,"chr"]==my_chr & mm10_known_circRNA[,"start"]>=my_start & mm10_known_circRNA[,"end"]<=my_end),]
  sub_merged <- mm10_merged_known_circRNA[(mm10_merged_known_circRNA[,"chr"]==my_chr & mm10_merged_known_circRNA[,"start"]>=my_start & mm10_merged_known_circRNA[,"end"]<=my_end),]
  
  
  result[i,4] <- my_end - my_start
  result[i,5] <- nrow(sub_Rybak)
  result[i,6] <- nrow(sub_other)
  result[i,7] <- nrow(sub_merged)
}
colnames(result)[c(4,5,6,7)] <- c("length", "Rybak_circRNAs", "other_circRNAs", "common_circRNAs")

qplot(result$Rybak_circRNAs,
      geom="histogram",
      fill=I("red"), 
      col=I("red"),
      alpha=I(.2),
      binwidth=1,
      main="Distribution of known circRNAs inside the predicted circRNA location",
      xlab="Known circRNAs (CircBase-Rybak)",
      ylab="Predicted circRNAs using circExplorer2")
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/circRNA/circRNA_prediction_vs_known_Rybak.png", sep = ""), width = 8, height = 4)


qplot(result$other_circRNAs,
      geom="histogram",
      fill=I("blue"), 
      col=I("blue"),
      alpha=I(.2),
      binwidth=1,
      main="Distribution of known circRNAs inside the predicted circRNA location",
      xlab="Known circRNAs (CircBase w/o Rybak)",
      ylab="Predicted circRNAs using circExplorer2")
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/circRNA/circRNA_prediction_vs_known_other.png", sep = ""), width = 8, height = 4)


qplot(result$common_circRNAs,
      geom="histogram",
      fill=I("blue"), 
      col=I("blue"),
      alpha=I(.2),
      binwidth=1,
      main="Distribution of known circRNAs inside the predicted circRNA location",
      xlab="Known circRNAs (common in both samples)",
      ylab="Predicted circRNAs using circExplorer2")
ggsave(paste("/data/home/students/ciora/circRNA-detection/plots/circRNA/circRNA_prediction_vs_common.png", sep = ""), width = 8, height = 4)


mm10_known_circRNA$other <- 1
mm10_Rybak_circRNA$Rybak <- 1
mm10_all_known_circRNA <- merge(mm10_known_circRNA, mm10_Rybak_circRNA, by = c("chr", "start", "end"), all = T)
is.na(mm10_all_known_circRNA) <- 0

png("/data/home/students/ciora/circRNA-detection/plots/circRNA/Venn_known_circRNAs_Rybak_vs_other.png", width = 400, height = 400)
grid.newpage()
draw.pairwise.venn(nrow(subset(mm10_all_known_circRNA, other == 1)), 
                   nrow(subset(mm10_all_known_circRNA, Rybak == 1)),
                   nrow(subset(mm10_all_known_circRNA, other == 1 & Rybak == 1)),
                   scaled = F,
                   category = c("other known circRNA", "Rybak circRNA"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
dev.off()


