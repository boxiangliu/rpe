#!/usr/bin/env Rscript
# bosh liu
# durga
# combine gender, PEER factors and genotype PCs into one covariate file

# library:
library('R.utils')
library('XLConnect')
library(data.table)
library(dplyr)
library(stringr)

# command arguments: 
args=commandArgs(trailingOnly=T,asValues=T,defaults=list(num_geno_pc=3,num_peer_factor=15))
genotype_pc_file=args$genotype_pc
peer_file=args$peer
gender_file=args$gender
output_file=args$output
gender_coding=args$gender_coding
num_geno_pc=args$num_geno_pc
num_peer_factor=args$num_peer_factor
row_and_colnames=as.logical(args$row_and_colnames)



# genotype_pc_file='../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv'
# peer_file='../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.PCs'
# gender_file='../processed_data/sex/sex/gender.tsv'
# output_file='../processed_data/sqtl/optimal_covariate/test_covariate/test.tsv'




# read input: 
message('INFO - reading input...')
genotype_pc=t(read.table(genotype_pc_file,header=T,row.names=1,check.names=F))
peer=read.table(peer_file,header=T,row.names=1,check.names=F)
gender=fread(gender_file,colClasses=c('character','character'))

dna2rna_fn='../data/meta/dna2rna.txt'
dna2rna=fread(dna2rna_fn,colClasses=c('character','character'))
gender=merge(gender,dna2rna,by.x='sample',by.y='DNA')
gender[,c('sample','RNA'):=list(RNA,NULL)]


# recode gender (M=0, F=1):
if (gender_coding=='numerical'){
	gender[,gender:=ifelse(gender=="M",0,1)]
} else if (gender_coding=='letter'){
	# do nothing
} else {stop("gender coding should be either 'numerical' or 'letter'")}


# Converting gender to data.frame:
gender=t(gender)
colnames(gender)=gender['sample',]
gender=as.data.frame(gender)
gender=gender[2,]


# check whether columns are in the same order:
message('INFO - checking column orders...')

stopifnot(colnames(peer)%in%colnames(genotype_pc))
if(!all(colnames(genotype_pc)==colnames(peer))){
	warning('Columns in genotype PC and PEER in different order.\nChanging the order of genotype PC to match PEER.')
	genotype_pc=genotype_pc[,match(colnames(peer),colnames(genotype_pc))]
}

stopifnot(colnames(peer)%in%colnames(gender))
if(!all(colnames(peer)==names(gender))){
	warning('Columns in gender and PEER in different order.\nChanging the order of gender to match PEER.')
	gender=gender[,match(colnames(peer),colnames(gender))]
}

stopifnot(all.equal(colnames(genotype_pc),colnames(peer)),all.equal(colnames(peer),names(gender)))


# combine covariates:
message('INFO - combining covariates...')
covariates=rbind(genotype_pc[1:num_geno_pc,],peer[1:num_peer_factor,],gender)


# set column names:
rownames(covariates)=c(paste0("C",seq(num_geno_pc)),paste0('InferredCov',seq(num_peer_factor)),'gender')
covariates=as.data.table(covariates,keep.rownames=T)
setnames(covariates,'rn','id')
# covariates[,ID:=NULL]


# transpose covariates:
# covariates=t(covariates)


# write output:
message('INFO - checking output...')
if (!dir.exists(dirname(output_file))){dir.create(dirname(output_file),recursive=TRUE)}
if (row_and_colnames){
	fwrite(covariates,file=output_file,sep='\t')
} else {
	fwrite(covariates,file=output_file,sep='\t')
}
