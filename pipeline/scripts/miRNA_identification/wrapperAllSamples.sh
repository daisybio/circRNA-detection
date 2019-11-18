#!/bin/bash
dataset=~/testing_data/small_testing/dataset_whole.txt
adapter="TGGAATTCTCGGGTGCCAAGG"
species=mmu
ref_dir=~/testing_data/reference/
ref_prefix=~/testing_data/reference/mm10
out_dir=~/testing_data/small_testing/output_whole/

bash ~/circRNA-detection/pipeline/scripts/miRNA_identification/miRNAmappingAllSamples.sh $dataset $adapter $species $ref_dir $ref_prefix $out_dir