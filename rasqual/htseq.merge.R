# Library: 
library(data.table)
library(stringr)


# Function: 
read_and_merge=function(filenames,suffix,id){
	y=data.table()
	for (fname in filenames){
		sample=basename(fname)
		sample=str_replace(sample,suffix,'')
		x=fread(fname)
		setnames(x,c(id,sample))
		if (nrow(y)!=0){
			y=merge(y,x,by=id)
		} else {
			y=x
		}
	}
	return(y)
}

# Read file names: 
glu_gene_files=list.files(path='../data/rnaseq/count/glucose/',pattern='*.glucose.gene_count',full.names=T)
glu_exon_files=list.files(path='../data/rnaseq/count/glucose/',pattern='*.glucose.exon_count',full.names=T)
gla_gene_files=list.files(path='../data/rnaseq/count/galactose/',pattern='*.galactose.gene_count',full.names=T)
gla_exon_files=list.files(path='../data/rnaseq/count/galactose/',pattern='*.galactose.exon_count',full.names=T)


# Read and merge files: 
glu_gene=read_and_merge(glu_gene_files,'.gene_count','gene_id')
glu_exon=read_and_merge(glu_exon_files,'.exon_count','exon_id')
gla_gene=read_and_merge(gla_gene_files,'.gene_count','gene_id')
gla_exon=read_and_merge(gla_exon_files,'.exon_count','exon_id')


# Remove special row: 
glu_gene=glu_gene[!str_detect(gene_id,'__'),]
glu_exon=glu_exon[!str_detect(exon_id,'__'),]
gla_gene=gla_gene[!str_detect(gene_id,'__'),]
gla_exon=gla_exon[!str_detect(exon_id,'__'),]


# Merge: 
gene=merge(glu_gene,gla_gene,by='gene_id')
exon=merge(glu_exon,gla_exon,by='exon_id')


# Write: 
if (!dir.exists('../data/rnaseq/count/merged/')) {dir.create('../data/rnaseq/count/merged/')}
fwrite(gene,'../data/rnaseq/count/merged/rpe.gene_count',sep='\t')
fwrite(exon,'../data/rnaseq/count/merged/rpe.exon_count',sep='\t')

