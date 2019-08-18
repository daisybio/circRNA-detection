#!/bin/bash
files=~/data/mouse_brain_GSE100265/fastq/miRNA_SRP096019/*
out_dir=/data/home/students/ciora/methods/miRDeep2/output/
ref_dir=/data/home/students/ciora/methods/miRDeep2/reference/
progress_list="/data/home/students/ciora/methods/miRDeep2/output/identification_progress.txt"
first=$1
last=$2
counter=0
processed=0
cd $out_dir
echo "Starting preprocessing"
echo "$(date) Starting preprocessing" >> $progress_list
for f in $files
do
	counter=$(($counter+1))
	if [ $counter -le $last ] && [ $first -le $counter ]; then
		sample=$(basename $f .fastq)
		echo "Started sample $counter $sample"
		echo "$(date) Started sample $counter $sample" >> $progress_list	
		reads_collapsed="${out_dir}mapping/${sample}_reads_collapsed.fa"
		arf="${out_dir}mapping/${sample}_reads_vs_ref.arf"
		mkdir "${out_dir}identification/${sample}/"
		cd "${out_dir}identification/${sample}/"
		miRDeep2.pl $reads_collapsed "${ref_dir}mm10.fa" $arf "${ref_dir}mature_ref.fa" "${ref_dir}mature_other.fa" "${ref_dir}hairpin_ref.fa" -t mmu 2>"${out_dir}identification/${sample}/${sample}_report.log"
		echo "Done sample $counter $sample"
		echo "$(date) Done sample $counter $sample" >> $progress_list
		processed=$(($processed+1))
	fi
done
echo "DONE RUN from $first to $last and processed $processed samples"
echo "$(date) DONE RUN from $first to $last and processed $processed samples" >> $progress_list
