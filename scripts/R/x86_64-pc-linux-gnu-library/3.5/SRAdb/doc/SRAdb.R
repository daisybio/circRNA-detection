### R code from vignette source 'SRAdb.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(width=50)


###################################################
### code chunk number 2: SRAdb.Rnw:56-58
###################################################
library(SRAdb)
sqlfile <- file.path(system.file('extdata', package='SRAdb'), 'SRAmetadb_demo.sqlite')


###################################################
### code chunk number 3: SRAdb.Rnw:63-68 (eval = FALSE)
###################################################
## 
## 
## timeStart <- proc.time()
## sqlfile <- getSRAdbFile()
## proc.time() - timeStart


###################################################
### code chunk number 4: SRAdb.Rnw:73-74
###################################################
file.info(sqlfile)


###################################################
### code chunk number 5: SRAdb.Rnw:79-80
###################################################
sra_con <- dbConnect(SQLite(),sqlfile)


###################################################
### code chunk number 6: SRAdb.Rnw:90-92
###################################################
sra_tables <- dbListTables(sra_con)
sra_tables


###################################################
### code chunk number 7: SRAdb.Rnw:96-97
###################################################
dbListFields(sra_con,"study")


###################################################
### code chunk number 8: SRAdb.Rnw:101-102
###################################################
dbGetQuery(sra_con,'PRAGMA TABLE_INFO(study)')


###################################################
### code chunk number 9: SRAdb.Rnw:106-108
###################################################
colDesc <- colDescriptions(sra_con=sra_con)[1:5,]
colDesc[, 1:4]


###################################################
### code chunk number 10: j1
###################################################
rs <- dbGetQuery(sra_con,"select * from study limit 3")
rs[, 1:3]


###################################################
### code chunk number 11: j2
###################################################
rs <- dbGetQuery(sra_con, paste( "select study_accession, 
        study_title from study where",
       "study_description like 'Transcriptome%'",sep=" "))
rs[1:3,]


###################################################
### code chunk number 12: SRAdb.Rnw:129-135
###################################################
getTableCounts <- function(tableName,conn) {
  sql <- sprintf("select count(*) from %s",tableName)
  return(dbGetQuery(conn,sql)[1,1])
}
do.call(rbind,sapply(sra_tables[c(2,4,5,11,12)], 
	getTableCounts, sra_con, simplify=FALSE))


###################################################
### code chunk number 13: SRAdb.Rnw:140-144
###################################################
rs <- dbGetQuery(sra_con, paste( "SELECT study_type AS StudyType, 
	count( * ) AS Number FROM `study` GROUP BY study_type order 
	by Number DESC ", sep=""))
rs


###################################################
### code chunk number 14: SRAdb.Rnw:149-153
###################################################
rs <- dbGetQuery(sra_con, paste( "SELECT instrument_model AS
	'Instrument Model', count( * ) AS Experiments FROM `experiment`
	GROUP BY instrument_model order by Experiments DESC", sep=""))
rs


###################################################
### code chunk number 15: SRAdb.Rnw:157-161
###################################################
rs <- dbGetQuery(sra_con, paste( "SELECT library_strategy AS 
	'Library Strategy', count( * ) AS Runs FROM `experiment` 
	GROUP BY library_strategy order by Runs DESC", sep=""))
rs


###################################################
### code chunk number 16: SRAdb.Rnw:169-171
###################################################
conversion <- sraConvert( c('SRP001007','SRP000931'), sra_con = sra_con )
conversion[1:3,]


###################################################
### code chunk number 17: SRAdb.Rnw:175-176
###################################################
apply(conversion, 2, unique)


###################################################
### code chunk number 18: SRAdb.Rnw:184-197
###################################################
rs <- getSRA( search_terms = "breast cancer", 
	out_types = c('run','study'), sra_con )
dim(rs)

rs <- getSRA( search_terms = "breast cancer", 
	out_types = c("submission", "study", "sample", 
	"experiment", "run"), sra_con )

# get counts for some information interested
apply( rs[, c('run','sample','study_type','platform',
	'instrument_model')], 2, function(x) 
	{length(unique(x))} )



###################################################
### code chunk number 19: SRAdb.Rnw:201-204
###################################################
rs <- getSRA (search_terms ='"breast cancer"',
	out_types=c('run','study'), sra_con)
dim(rs)


###################################################
### code chunk number 20: SRAdb.Rnw:208-211
###################################################
rs <- getSRA( search_terms ='MCF7 OR "MCF-7"',
	out_types = c('sample'), sra_con ) 
dim(rs)


###################################################
### code chunk number 21: SRAdb.Rnw:215-218
###################################################
rs <- getSRA( search_terms ='submission_center: GEO', 
     out_types = c('submission'), sra_con )  
dim(rs)


###################################################
### code chunk number 22: SRAdb.Rnw:222-225
###################################################
rs <- getSRA( search_terms ='Carcino*', 
     out_types = c('study'), sra_con=sra_con )  
dim(rs)


###################################################
### code chunk number 23: SRAdb.Rnw:232-233
###################################################
rs = listSRAfile( c("SRX000122"), sra_con, fileType = 'sra' )


###################################################
### code chunk number 24: SRAdb.Rnw:238-240
###################################################
rs = getSRAinfo ( c("SRX000122"), sra_con, sraType = "sra" )
rs[1:3,]


###################################################
### code chunk number 25: SRAdb.Rnw:245-246 (eval = FALSE)
###################################################
## getSRAfile( c("SRR000648","SRR000657"), sra_con, fileType = 'sra' )


###################################################
### code chunk number 26: SRAdb.Rnw:252-253 (eval = FALSE)
###################################################
## system ("fastq-dump SRR000648.sra")


###################################################
### code chunk number 27: SRAdb.Rnw:257-260 (eval = FALSE)
###################################################
## getFASTQinfo( c("SRR000648","SRR000657"), sra_con, srcType = 'ftp' )
## 
## getSRAfile( c("SRR000648","SRR000657"), sra_con, fileType = 'fastq' )


###################################################
### code chunk number 28: SRAdb.Rnw:267-287 (eval = FALSE)
###################################################
## ## List fasp addresses for associated fastq files:
## listSRAfile ( c("SRX000122"), sra_con, fileType = 'fastq', srcType='fasp')
## 
## ## get fasp addresses for associated fastq files:
## getFASTQinfo( c("SRX000122"), sra_con, srcType = 'fasp' )
## 
## ## download fastq files using fasp protocol:
## # the following ascpCMD needs to be constructed according custom 
## # system configuration
## # common ascp installation in a Linux system:
## ascpCMD <-  'ascp -QT -l 300m -i 
##  /usr/local/aspera/connect/etc/asperaweb_id_dsa.putty'
## 
## ## common ascpCMD for a Mac OS X system:
## # ascpCMD <- "'/Applications/Aspera Connect.app/Contents/
## #  Resources/ascp' -QT -l 300m -i '/Applications/
## # Aspera Connect.app/Contents/Resources/asperaweb_id_dsa.putty'"
##    
## getSRAfile( c("SRX000122"), sra_con, fileType = 'fastq', 
## 	srcType = 'fasp',  ascpCMD = ascpCMD )


###################################################
### code chunk number 29: SRAdb.Rnw:291-297 (eval = FALSE)
###################################################
## ## List fasp addresses of sra files associated with "SRX000122"
## listSRAfile( c("SRX000122"), sra_con, fileType = 'sra', srcType='fasp')
## 
## ## download sra files using fasp protocol
## getSRAfile( c("SRX000122"), sra_con, fileType = 'sra', 
##  srcType = 'fasp',  ascpCMD = ascpCMD )


###################################################
### code chunk number 30: SRAdb.Rnw:312-313 (eval = FALSE)
###################################################
## startIGV("mm")


###################################################
### code chunk number 31: SRAdb.Rnw:317-324 (eval = FALSE)
###################################################
## exampleBams = file.path(system.file('extdata',package='SRAdb'),
##   dir(system.file('extdata',package='SRAdb'),pattern='bam$'))
## sock <- IGVsocket()
## IGVgenome(sock, 'hg18')
## IGVload(sock, exampleBams)
## IGVgoto(sock, 'chr1:1-1000')
## IGVsnapshot(sock)


###################################################
### code chunk number 32: SRAdb.Rnw:339-351 (eval = FALSE)
###################################################
## library(SRAdb)
## library(Rgraphviz)
## 
## g <- sraGraph('primary thyroid cell line', sra_con)
## attrs <- getDefaultAttrs(list(node=list(
## 	fillcolor='lightblue', shape='ellipse')))
## plot(g, attrs=attrs)
## 
## ## similiar search as the above, returned much larger data.frame and graph is too clouded
## g <- sraGraph('Ewing Sarcoma', sra_con)
## plot(g)	
## 


###################################################
### code chunk number 33: SRAdb.Rnw:358-359
###################################################
dbDisconnect(sra_con)


###################################################
### code chunk number 34: SRAdb.Rnw:371-390 (eval = FALSE)
###################################################
## 
## library(SRAdb)
## 
## setwd('1000g')
## if( ! file.exists('SRAmetadb.sqlite') ) {
## 	sqlfile <- getSRAdbFile()
## } else {
## 	sqlfile <- 'SRAmetadb.sqlite'	
## }
## sra_con <- dbConnect(SQLite(),sqlfile)
## 
## ## get all related accessions
## rs <- getSRA( search_terms = '"1000 Genomes Project"', 
## 	sra_con=sra_con, acc_only=TRUE)
## dim(rs)
## head(rs)
## 
## ## get counts for each data types
## apply( rs, 2, function(x) {length(unique(x))} )


###################################################
### code chunk number 35: SRAdb.Rnw:395-397 (eval = FALSE)
###################################################
## runs <- tail(rs$run)
## fs <- getSRAinfo( runs, sra_con, sraType = "sra" )


###################################################
### code chunk number 36: SRAdb.Rnw:401-402 (eval = FALSE)
###################################################
## getSRAfile( runs, sra_con, fileType ='sra', srcType = "ftp" )


###################################################
### code chunk number 37: SRAdb.Rnw:406-409 (eval = FALSE)
###################################################
## ascpCMD <- "'/Applications/Aspera Connect.app/Contents/Resources/ascp' -QT -l 300m -i '/Applications/Aspera Connect.app/Contents/Resources/asperaweb_id_dsa.putty'"
## 
## sra_files = getSRAfile( runs, sra_con, fileType ='sra', srcType = "fasp", ascpCMD = ascpCMD )


###################################################
### code chunk number 38: SRAdb.Rnw:413-416 (eval = FALSE)
###################################################
## for( fq in basename(sra_files$fasp) ) {
## 	system ("fastq-dump SRR000648.lite.sra")
## }


###################################################
### code chunk number 39: SRAdb.Rnw:424-425
###################################################
toLatex(sessionInfo())


