#!/bin/bash
files=~/data/mouse_brain_GSE100265/fastq/miRNA_SRP096019/*
cd ~/data/mouse_brain_GSE100265/fastq/miRNA_clipped
for f in $files
do
	sample=$(basename $f .fastq)
	flexbar -r $f -a ~/circRNA-detection/scripts/adapters_miRNA.fa --adapter-relaxed --adapter-trim-end RIGHT --adapter-gap -3 --min-read-length 17 -t $sample
done






