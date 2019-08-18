#/bin/bash
sraFiles=~/data/mouse_brain_GSE100265/raw/mRNA_SRP096017/*.sra
for f in $sraFiles
do
	  echo "Processing $f file..."
	~/circRNA-detection/soft/sratoolkit.2.9.6-1-centos_linux64/bin/fastq-dump -I --split-files $f
done
