library(data.table)
library(cowplot)

GEDi_fn = '../data/eye_disease/GEDi.txt'
zscore_fn = '../processed_data/rpe_specific_genes.GTExV7/all_genes.txt'

read_GEDi = function(GEDi_fn){
	GEDi = fread(GEDi_fn,header=TRUE)
	return(GEDi)
}


read_zscore = function(zscore_fn){
	zscore = fread(zscore_fn)
	return(zscore)
}


GEDi = read_GEDi(GEDi_fn)
zscore = read_zscore(zscore_fn)
zscore = zscore[type=='protein_coding'&chr%in%paste0('chr',1:22)]

zscore[,GEDi := gene_name %in% GEDi$gene_name]

ggplot(zscore,aes(GEDi,zscore)) + geom_violin()
ggplot(zscore,aes(GEDi,zscore)) + geom_boxplot()
ggplot(zscore,aes(zscore,color=GEDi)) + geom_density()

=
