#!/bin/bash
file_list=$1
output_dir=$2
index_dir=$3
if [ -d "$output_dir" ];then
	echo "$output_dir exists"
#	rm -r "$output_dir"
else
	mkdir "$output_dir"
fi
progress_file="$output_dir/processed_files.txt"
#touch $progress_file
while read line
do
	sample=$(basename $line .fastq)
	echo "Processing $sample"
	echo "Processing $sample" >> $progress_file
	sample_out="$output_dir/$sample/"
	mkdir "$sample_out"
	fastq=$line
	echo "Starting Salmon"
	echo "$(date) $sample Starting Salmon" >> $progress_file
	salmon quant -i $index_dir --libType A -r $fastq --validateMappings -o $sample_out
#	salmon quant -i $index_dir --libType SR -r $fastq --validateMappings -o $sample_out
	echo "$sample done"
	echo "$(date) $sample done" >> $progress_file
done < $file_list
echo "DONE!!!"
