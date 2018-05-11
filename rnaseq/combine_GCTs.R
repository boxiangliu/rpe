library(data.table)
library(foreach)
library(stringr)

args = commandArgs(T)
in_fn = args[1]
out_fn = args[2]

file_list = fread(in_fn,header=FALSE)
counts = foreach(i = 1:nrow(file_list),.combine = function(x,y) merge(x,y,all=TRUE,by=c('Name','Description'))) %do% {
	fn = file_list[i,V1]
	print(fn)
	x = fread(fn,header=TRUE,skip=2)
	return(x)
}
fwrite(counts,out_fn,sep='\t')