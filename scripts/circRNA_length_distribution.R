circRNAs_unfiltered <- read.table("/nfs/home/students/ciora/methods/circexplorer2_method/circRNA_results.tsv", sep="\t", header=T, stringsAsFactors = F)
circRNAs_unfiltered <- circRNAs_unfiltered[,c(1,2,3,4)]

circRNAs_unfiltered$length <- circRNAs_unfiltered$stop - circRNAs_unfiltered$start


library(scales)

ggplot(circRNAs_unfiltered, aes(x=length)) +
  geom_histogram(colour=I("red"), fill=I("red"), alpha=I(.2)) +
  labs(title = "circRNA length distribution (unfiltered)") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

ggsave("/nfs/home/students/ciora/plots/analysis/circRNA_length_distribution_unfiltered.png", width = 6, height = 4)


f = function(x, output) {
  # x is the row of type Character
  circRNA = as.character(x[1])
  start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA),':'), "[", 2),'-'), "[", 1))
  end <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(circRNA),'-'), "[", 2),'_'), "[", 1))
  length <- end - start
  freq = x[2]
  print(length)
}


labeled_circRNAs = read.table("/nfs/home/students/ciora/methods/results_analysis/circRNA_ATXN1_notation.txt", header = T, stringsAsFactors = F)
labeled_circRNAs$length <- labeled_circRNAs$end - labeled_circRNAs$start
labeled_circRNAs$y_position <- c(60,60,60,60,60,60,60, 60)

ggplot(circRNAs_unfiltered, aes(x=length)) +
  geom_histogram(colour=I("red"), fill=I("red"), alpha=I(.2)) +
  labs(title = "circRNA length distribution + ATXN1 circRNAs") + 
  geom_vline(xintercept = labeled_circRNAs$length, color = "black", size=0.7) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data = labeled_circRNAs, mapping = aes(label = notation, y = y_position, x = length), angle = 90, vjust = - 0.3)
ggsave("/nfs/home/students/ciora/plots/analysis/circRNA_length_distribution_and_ATXN1.png", width = 6, height = 4)
