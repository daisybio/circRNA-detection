# circRNA Sponging Analysis on Mouse Data

## Overview
This analysis consists of several setps:
1. Detection of circRNAs from total RNA-seq data
2. Quantification of miRNAs from small RNA-seq data
3. Quantification of mRNAs from RNA-ses data
4. Detection of miRNA binding sites on previously found circRNAs
5. Correlation analysis between miRNAs and circRNAs
6. Further analysis

## Dataset
This analysis was done on mouse data covering four brain regions (cerebellum, cortex, hippocampus, olfactory bulb) and two conditions (wild-type and knockout). This data was obtained through three different sequencing methods (RNA-seq, total RNA-seq, small RNA-seq) The dataset was produced by the Rajewsky lab in 2017 and is available in the GEO database under the accession number GSE93130.

## Documentation
The ```documentation/``` folder contains step-by-step descriptions of the quantification of miRNAs, circRNAs, mRNAs and detection of binding sites. These are provided as Jupyter Notebooks and include the installation and usage of external tools.

## Results
The ```results/``` folder contains the expression data for miRNAs, mRNAs, circRNAs and correlation results.

## Scripts
All used scripts are provided in the ```scripts/``` folder.
```bash
addFilesToList.sh                       # get a list of all files in a folder
amount_circRNA.R                        # visualization of circRNA detection results        
amount_mRNA.R                           # visualization of circRNA detection results                                   
binding_sites_plot.R                    # analysis and visualization of binding sites results           
circRNA_length_distribution.R           # visualization of lengths for detected circRNAs         
clipAdapters.sh                         # clip adapters using flexbar
correlation.R                           # correlation analysis for all circRNA-miRNA pairs and vizualisation of results
countAdapters.sh                        # check whether adapter clipping was successfull
DiffExp_miRNA.R                         # differential expression analysis for miRNA data (WT vs. KO)
downloadSRA.R                           # download all samples of a dataset in SRA format
genome_vizualisation.R                  # genomic visualization of sponged miRNAs and circRNAs found on the ATXN1 gene
get_circRNA_sequences.sh                # extract fastq sequences for all detected circRNAs
miRNA_comparison.R                      # consistency analysis for miRNAs
miRNA_identification_miRDeep2.sh        # miRNA quantification for a list of samples unsing miRDeep2
mRNA_to_gene_level.R                    # transform transcript-level counts to gene-level counts for mRNA Salmon mapping
preprocesing_miRNA_miRDeep2.sh          # data preparation for miRDeep2
runCIRCExplorer2ForList.sh              # circRNA detection for a list of samples using CIRCExplorer2
runSalmonForList.sh                     # mRNA quantification for a list of samples using Salmon
samples_PCA.R                           # PCA sample clustering and visualization of mRNA counts
SRAtoFastq.sh                           # convert SRA files to fastq
```



