#!/bin/bash
file_list=$1
output_dir=$2
reference_dir=$3
#if [ -d "$output_dir" ];then
#	rm -r "$output_dir"
#fi
#mkdir "$output_dir"
progress_file="$output_dir/processed_files.txt"
#touch $progress_file
while read line
do
	sample=$(basename $line .fastq)
	echo "Processing $sample"
	echo "Processing $sample" >> $progress_file
	sample_out="$output_dir/$sample/"
	mkdir "$sample_out"
	sample_out="$sample_out/circ_out"
	mkdir "$sample_out"
	fastq=$line
	echo "Starting TopHat2"
	echo "$(date) $sample Starting TopHat2" >> $progress_file
	tophat2 -a 6 --microexon-search -m 2 -p 16 -G "$reference_dir/mm10_kg.gtf" -o "$sample_out/tophat/" "$reference_dir/mm10" $line
	echo "BamToFastq"
	echo "$(date) $sample BamToFastq" >> $progress_file
	bamToFastq -i "$sample_out/tophat/unmapped.bam" -fq "$sample_out/tophat/unmapped.fastq"
	echo "Starting TopHat Fusion"
	echo "$(date) $sample Starting TopHat Fusion" >> $progress_file
	tophat2 -o "$sample_out/tophat_fusion" -p 15 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search "$reference_dir/mm10" "$sample_out/tophat/unmapped.fastq"
	echo "Starting CIRCexplorer2 parse"
	echo "$(date) $sample Starting CIRCexplorer2 parse" >> $progress_file
	CIRCexplorer2 parse -b $sample_out/back_spliced_junction.bed -t TopHat-Fusion "$sample_out/tophat_fusion/accepted_hits.bam" > CIRCexplorer2_parse.log
	echo "Starting CIRCexplorer2 annotate"
	echo "$(date) $sample Starting CIRCexplorer2 annotate" >> $progress_file
	CIRCexplorer2 annotate -r "$reference_dir/mm10_ref.txt" -g "$reference_dir/mm10.fa" -b $sample_out/back_spliced_junction.bed -o $sample_out/circularRNA_known.txt > CIRCexplorer2_annotate.log
	echo "$(date) $sample Annotation done" >> $progress_file
	echo "Annotation done"
	echo "$sample done" >> $progress_file
done < $file_list
echo "DONE!!!"
