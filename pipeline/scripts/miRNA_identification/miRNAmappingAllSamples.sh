#!/bin/bash
dataset=$1
adapter=$2
species=$3
ref_dir=$4
ref_prefix=$5
out_dir=$6

mkdir -p $out_dir

while read line
do
	sample=$(cut -f1 <<< $line)
	miRNAfastq=$(cut -f3 <<< $line)
	bash ~/circRNA-detection/pipeline/scripts/miRNAidentification/miRNAmappingForSample.sh $sample $miRNAfastq $adapter $species $ref_dir $ref_prefix $out_dir
done < $dataset
