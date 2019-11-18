#!/bin/bash
dataset=~/testing_data/big_testing/dataset.tsv
adapter="TGGAATTCTCGGGTGCCAAGG"
species=mmu
ref_dir=~/testing_data/reference/
ref_prefix=~/testing_data/reference/mm10
out_dir=~/testing_data/big_testing/after_mapping/
scripts_dir=~/circRNA-detection/pipeline/scripts/

echo "Started pipeline"
mkdir -p "${out_dir}/samples/"

#run circRNA detection module
echo "circRNA detection module using CIRCexplorer2"
#bash ${scripts_dir}/circRNA_detection/circDetectionAllSamples.sh $dataset $ref_prefix ${out_dir}/samples
mkdir -p "${out_dir}/results/circRNA"
Rscript ${scripts_dir}/circRNA_detection/circRNA_results_processing.R $dataset "${out_dir}"

#run miRNA identification module
echo "circRNA identification module using miRDeep2"
#bash ${scripts_dir}/miRNA_identification/miRNAmappingAllSamples.sh $dataset $adapter $species $ref_dir $ref_prefix $out_dir
mkdir -p "${out_dir}/results/miRNA/read_distributions"
Rscript ${scripts_dir}/miRNA_identification/miRNA_results_processing.R $dataset "${out_dir}"

echo "Extracting fasta sequences for circRNAs"
mkdir -p "${out_dir}/results/binding_sites/input/"
bash "${scripts_dir}/binding_sites/get_circRNA_sequences.sh" $ref_prefix $out_dir

echo "Detection of miRNA binding sites on circRNAs using miRanda"
mkdir -p "${out_dir}/results/binding_sites/output/"
~/methods/miRanda/miRanda_aug_2010/bin/miranda "${ref_dir}/mature_ref.fa" "${out_dir}/results/binding_sites/input/circRNAs.fa" -out "${out_dir}/results/binding_sites/output/bind_sites_raw.out" -quiet

echo "Processing binding sites"
echo -e "miRNA\tTarget\tScore\tEnergy-Kcal/Mol\tQuery-Aln(start-end)\tSubject-Al(Start-End)\tAl-Len\tSubject-Identity\tQuery-Identity" > "${out_dir}/results/binding_sites/output/circRNA_bind_sites_results.txt"
grep -A 1 "Scores for this hit:" "${out_dir}/results/binding_sites/output/bind_sites_raw.out" | sort | grep ">" | cut -c 2- >> "${out_dir}/results/binding_sites/output/circRNA_bind_sites_results.txt"

mkdir -p "${out_dir}/results/binding_sites/plots/"
echo "Filter and analyze binding sites"
Rscript "${scripts_dir}/binding_sites/binding_sites_analysis.R" $out_dir

mkdir -p "${out_dir}/results/sponging/plots/"
echo "Sponging analysis"
Rscript "${scripts_dir}/sponging/correlation.R" $dataset $out_dir

echo "DONE!"

