library(data.table)
library(foreach)
library(stringr)

args = commandArgs(T)
in_fn = args[1]
out_fn = args[2]

file_list = fread(in_fn,header=FALSE)
counts = foreach(i = 1:nrow(file_list),.combine = function(x,y) merge(x,y,all=TRUE,by=c('Transcript','Gene_Name'))) %do% {
	fn = file_list[i,V1]
	sample = str_replace(basename(fn),'.metrics.tmp.txt.intronReport.txt','')
	print(sample)
	x = fread(fn,header=TRUE)[,1:3]
	setnames(x,'Exon_Reads',sample)
	return(x)
}
counts[is.na(counts)]=0
setnames(counts,c('Transcript','Gene_Name'),c('Name','Description'))
fwrite(counts,out_fn,sep='\t')