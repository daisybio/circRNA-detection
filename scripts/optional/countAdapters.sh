#!/bin/bash
adapter="TGGAATTC"
files=~/data/mouse_brain_GSE100265/fastq/miRNA_SRP096019/*
list="adapter_counts.tsv"
echo -e "sample\t#adapters\t#sequences\tpercent" > $list
for f in $files
do
	sample=$(basename $f .fastq)
	adapters=$(grep -c $adapter $f)
	total=$(grep -c @SRR $f)
	if [ $total -eq 0 ]; then
		percent=0
	else
		percent=$((adapters*100/total))
	fi	
	echo -e "$sample\t$adapters\t$total\t$percent" >> $list
done






