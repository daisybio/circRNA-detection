#!/bin/bash
if [ -z $1 ]
then files=~/data/mouse_brain_GSE100265/fastq/rRNAdepleted_SRP109830/*
else files=$1/*
fi
if [ -z $2 ]
then list=~/circRNA-detection/scripts/fastq_list_rRNAdeplteted.txt
else list=$2
fi
rm $list
for f in $files
do
	echo $f >> $list
done
