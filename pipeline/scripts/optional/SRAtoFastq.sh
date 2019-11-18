#/bin/bash
sraFiles=~/testing_data/sra/totalRNA/*.sra
output_dir=~/testing_data/fastq/totalRNA/
for f in $sraFiles
do
	  echo "Processing $f file..."
	~/soft/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump -I --split-files -O $output_dir $f
done
