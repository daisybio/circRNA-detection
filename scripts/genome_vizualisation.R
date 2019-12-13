sample_correspondence <- read.table("/nfs/home/students/ciora/data/mouse_brain_GSE100265/circRNA_miRNA_mRNA.tsv", header = T, stringsAsFactors = F, sep = "\t")

raw_bindSites <- read.table("/nfs/home/students/ciora/circRNA-detection/results/miRanda/bindsites_25%_filtered.tsv", header = T, stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]
allBindSites <- count(bindSites, Target, name="freq")

zero_counts_sample_cutoff <- ceiling(0.2*nrow(dataset))

miRNA_expression_raw <- read.table("/nfs/home/students/ciora/circRNA-detection/results/miRNA/miRDeep2/miRNA_counts_all_samples.tsv", header = T, stringsAsFactors = F)
miRNA_expression <- miRNA_expression_raw[rowSums(miRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]

circRNA_expression_raw <- read.table("/nfs/home/students/ciora/circRNA-detection/results/circExplorer2/circRNA_filtered_results.tsv",header = T, stringsAsFactors = F)
circRNA_expression <- circRNA_expression_raw[rowSums(circRNA_expression_raw == 0) <= zero_counts_sample_cutoff, ]



