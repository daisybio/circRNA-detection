#source('http://bioconductor.org/biocLite.R')
#biocLite('SRAdb')
library(SRAdb)
srafile = getSRAdbFile()
file.info(srafile)
con = dbConnect(RSQLite::SQLite(),srafile)
#rRNA depleted
#getSRAfile('SRP109830',con,fileType='sra')
#mRNA
getSRAfile('SRP096017',con,fileType='sra')
#miRNA
getSRAfile('SRP096019',con,fileType='sra')