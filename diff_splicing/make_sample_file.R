# Make sample file
# Boxiang Liu
# 2017-12-11
library(data.table)
library(stringr)

glu_junc_fn='../data/rnaseq/leafcutter/glucose/juncfiles.txt'
gal_junc_fn='../data/rnaseq/leafcutter/galactose/juncfiles.txt'
geno_pc_fn='../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv'
gender_fn='../processed_data/sex/sex/gender.tsv'
out_dir='../processed_data/diff_splicing/sample_file/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_junc_file=function(fn){
	junc=fread(fn,header=FALSE,col.names=c('filename'))
	junc[,filename:=str_replace(filename,'.junc','')]
	junc[,sample:=basename(filename)]
	return(junc)
}

make_sample_table=function(x1,x2,group1,group2){
	x1=cbind(x1,treatment=group1)
	x2=cbind(x2,treatment=group2)
	return(rbind(x1,x2))
}


read_geno_pc=function(geno_pc_fn,num){
	geno_pc=fread(geno_pc_fn)
	geno_pc=geno_pc[,1:(num+1)]
	return(geno_pc)
}

read_gender=function(gender_fn){
	gender=fread(gender_fn,colClasses='character')

	dna2rna_fn='../data/meta/dna2rna.txt'
	dna2rna=fread(dna2rna_fn,colClasses=c('character','character'))
	gender=merge(gender,dna2rna,by.x='sample',by.y='DNA')
	gender[,c('sample','RNA'):=list(RNA,NULL)]

	return(gender)
}

add_covariate=function(sample_table,geno_pc,gender){
	if (!all(geno_pc$sample%in%sample_table$sample)){
		stop('genotype PC file and sample table have different sample names!')
	}

	if (!all(sample_table$sample%in%gender$sample)){
		stop('gender file and sample table have different sample names!')
	}

	merged=merge(sample_table,geno_pc,by='sample',all=FALSE)
	merged=merge(merged,gender,by='sample',all=FALSE)
	return(merged)
}

main=function(){
	glu_junc=read_junc_file(glu_junc_fn)
	gal_junc=read_junc_file(gal_junc_fn)
	sample_table=make_sample_table(glu_junc,gal_junc,'glu','gal')
	geno_pc=read_geno_pc(geno_pc_fn,2)
	gender=read_gender(gender_fn)
	sample_table=add_covariate(sample_table,geno_pc,gender)
	fwrite(sample_table[,list(filename,treatment,PC1,PC2,gender)],sprintf('%s/sample_table.txt',out_dir),sep='\t',col.names=FALSE)
}

main()