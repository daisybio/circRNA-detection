#!/bin/bash
sample=testing63
fastq=~/testing_data/short_fastq/miRNA/63_100reads.fastq
adapter="TGGAATTCTCGGGTGCCAAGG"
species=mmu
ref_dir=~/testing_data/reference/
ref_prefix=~/testing_data/reference/mm10
out_dir=~/testing_data/small_testing/output100/



bash ~/circRNA-detection/pipeline/scripts/miRNA_identification/miRNAmappingForSample.sh $sample $fastq $adapter $species $ref_dir $ref_prefix $out_dir