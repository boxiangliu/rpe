library(data.table)
library(dplyr)
library(dtplyr)
library(rasqualTools)

# Function: 
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


# Read HTSeq count matrix: 
count=fread('../data/rnaseq/count/merged/rpe.gene_count')
count_df=as.data.frame(count[,2:ncol(count)])
rownames(count_df)=unlist(count[,1])


# Subset to protein coding and lncRNA:
count_df=subset_to_protein_and_lncRNA(count_df)


# Remove unexpressed genes: 
count_df=remove_unexpressed(count_df)


# Split count matrix into glucose and galactose: 
glu=count_df%>%select(contains('lucose'))
gla=count_df%>%select(contains('lactose'))



# Save count matrix in both binary format and plain text format: 
saveRasqualMatrices(list(glucose=glu,galactose=gla),"../processed_data/rasqual/expression/", file_suffix = "expression")


# Read gene GC content: 
gc=fread('../processed_data/rasqual/gcc.exon.txt',header=F)
setnames(gc,c('gene_id','percentage_gc_content'))


# Calculate offset: 
size_factors_glu=rasqualCalculateSampleOffsets(glu,gc,gc_correct = TRUE)
saveRasqualMatrices(list(glucose=glu,galactose=gla), "../processed_data/rasqual/expression/", file_suffix = "size_factors_gc")


# Generate rasqual input file: 
for (i in 1:22){
	system(sprintf('grep chr%s ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf | python rasqual/make_input.py ../data/genotype/asvcf/glucose/rpe.imputed.chr%s.all_filters.vcf.new.gz 1000000 > ../processed_data/rasqual/input/rasqual.input.chr%s.txt',i,i,i))
}


# Filter rasqual input: 
for (i in 1:22){
	x=fread(sprintf('../processed_data/rasqual/input/rasqual.input.chr%s.txt',i))
	y=x[V1%in%rownames(count_df),]
	fwrite(y,sprintf('../processed_data/rasqual/input/rasqual.input.chr%s.filt.txt',i),sep='\t',col.names=F)
}
