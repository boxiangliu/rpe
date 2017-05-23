library(dplyr)
library(stringr)

mat=read.table('../data/rnaseq/rpkm/rpe.filt.quant_norm.rpkm',header=T,check.names=F)
rownames(mat)=mat[,1]
glu=mat%>%select(contains('glucose'))
glu$gene_id=rownames(glu)
gla=mat%>%select(contains('galactose'))
gla$gene_id=rownames(gla)

gene_annotation=read.table('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',header=F,stringsAsFactors=F,col.names=c('chr','start_bak','end_bak','strand','gene_id','type'))%>%mutate(chr=paste0('chr',chr),start=ifelse(strand=='+',start_bak,end_bak),end=ifelse(strand=='+',end_bak,start_bak))%>%select(chr,start,end,gene_id)

glu=merge(gene_annotation,glu,by='gene_id')
glu=glu%>%select(c(2:4,1,5:28))
glu=glu%>%arrange(chr,start)
colnames(glu)[1]='#chr'
colnames(glu)=str_replace(colnames(glu),'.glucose','')
write.table(glu,'../processed_data/fastqtl/expression/glucose.bed',quote=F,row.names=F,col.names=T,sep='\t')

gla=merge(gene_annotation,gla,by='gene_id')
gla=gla%>%select(c(2:4,1,5:28))
gla=gla%>%arrange(chr,start)
colnames(gla)[1]='#chr'
colnames(gla)=str_replace(colnames(gla),'.galactose','')
write.table(gla,'../processed_data/fastqtl/expression/galactose.bed',quote=F,row.names=F,col.names=T,sep='\t')