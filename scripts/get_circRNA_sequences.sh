#!/bin/bash
input=~/methods/circexplorer2_method/circRNA_filtered_results.tsv
#input=~/methods/circexplorer2_method/cdr1as.tsv
fasta=~/methods/circexplorer2_method/reference/mm10.fa
output=~/methods/miRanda/input/circRNA.fa
#output=~/methods/miRanda/input/cdr1as.test
rm $output
counter=0
while read line
do
	chr=$(echo $line | awk '{ print $1 }')
	from=$(echo $line | awk '{ print $2 }')
	to=$(echo $line | awk '{ print $3 }')
	strand=$(echo $line | awk '{ print $4 }')
	sam_out=$(samtools faidx $fasta "$chr:$from-$to")
	header=">${chr}:${from}-${to}_${strand}"
	echo $header >> $output
	if [ $strand == "+" ]; then
		sequence=$(echo "$sam_out" 2>&1 | tail -n +2 )
	elif [ $strand == "-" ]; then
		sequence=$(echo "$sam_out" 2>&1 | tail -n +2 | rev | tr {AGTCagtc} {TCAGtcag})

	fi		
		echo "$sequence" >> $output
	counter=$(($counter+1))
	if [ $(($counter%100)) -eq 0 ]; then
		echo $counter
	fi		

done < <(tail -n +2 $input)


