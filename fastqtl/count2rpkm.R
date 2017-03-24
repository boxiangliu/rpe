library(data.table)
library(dplyr)
# Functions: 
remove_unexpressed=function(x){
	y=rowSums(x>=6)
	return(x[y>=10,])
}

subset_to_protein_and_lncRNA=function(x){
	gencode=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed')[,5:6,with=F]
	setnames(gencode,c('gene_id','type'))
	y=x[rownames(x)%in%gencode[type=='lincRNA'|type=='protein_coding',gene_id],]
	return(y)
}

divide_by_gene_length=function(count,gene_length){
	gene_length=gene_length[match(rownames(count),gene_length$gene_id),]
	return((count/gene_length$length)*1000)
}

# Read: 
count=read.table('../data/rnaseq/count/merged/rpe.gene_count',header=T,check.names=F)
rownames(count)=count[,1]
count[,1]=NULL


# Calculate RPKM:
rpm=sweep(count,2,colSums(count)/1e6,'/')
gene_length=read.table('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',header=F,stringsAsFactors=F)%>%transmute(gene_id=V5,length=abs(V3-V2))
rpkm=divide_by_gene_length(rpm,gene_length)


# To output:
rpkm$gene_id=rownames(rpkm)
rpkm=rpkm[,c(ncol(rpkm),1:(ncol(rpkm)-1))]
write.table(rpkm,'../data/rnaseq/rpkm/rpe.rpkm',row.names=F,col.names=T,sep='\t',quote=F)

# Remove unexpressed genes: 
expressed_genes=rownames(remove_unexpressed(count))
rpkm=rpkm[rownames(rpkm)%in%expressed_genes,]

# Keep only protein coding and lncRNA:
rpkm=subset_to_protein_and_lncRNA(rpkm)

# To output:
write.table(rpkm,'../data/rnaseq/rpkm/rpe.filt.rpkm',row.names=F,col.names=T,sep='\t',quote=F)

