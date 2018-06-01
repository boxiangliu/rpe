library(data.table)
library(stringr)

gencode_fn = '../data/reference/gencode.v19.annotation.gtf'

read_gencode = function(gencode_fn){
	x = fread(gencode_fn,select=c(1,3,4,5,7,9))
	setnames(x,c('chr','type','start','end','strand','annotation'))
	x = x[type=='gene',]
	x[,gene_id:=str_extract(annotation,'(?<=gene_id \\")(.+?)(?=\\";)')]
	x[,gene_name:=str_extract(annotation,'(?<=gene_name \\")(.+?)(?=\\";)')]
	x[,gene_type:=str_extract(annotation,'(?<=gene_type \\")(.+?)(?=\\";)')]
	x$annotation=NULL
	return(x)
}

