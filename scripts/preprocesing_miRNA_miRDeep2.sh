#!/bin/bash
adapter="TGGAATTCTCGGGTGCCAAGG"
files=~/data/mouse_brain_GSE100265/fastq/miRNA_SRP096019/*
out_dir=/nfs/home/students/ciora/methods/miRDeep2/output/mapping/
ref_prefix=/nfs/home/students/ciora/methods/miRDeep2/reference/mm10
progress_list="/nfs/home/students/ciora/methods/miRDeep2/output/mapping_progress.txt"
echo "Starting preprocessing"
echo "Starting preprocessing" > $progress_list
for f in $files
do
	sample=$(basename $f .fastq)
	echo "Preprocessing sample $sample"
	echo "Preprocessing sample $sample" >> $progress_list	
	mapper.pl $f -e -h -i -j -k $adapter -l 18 -m -p $ref_prefix -s "${out_dir}${sample}_reads_collapsed.fa" -t "${out_dir}${sample}_reads_vs_ref.arf" -v -o 4 >> $progress_list
	echo "Preprocessing done for sample $sample"
	echo "Preprocessing done for sample $sample" >> $progress_list

done
