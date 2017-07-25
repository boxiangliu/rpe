library(data.table)
library(stringr)
library(dplyr)
library(dtplyr)

count_fn='../data/rnaseq/count/merged/rpe.gene_count'
glu_cov_fn='../processed_data/select_covariates/select_covariates/covariates/glu_cov.8.txt'
gal_cov_fn='../processed_data/select_covariates/select_covariates/covariates/gal_cov.9.txt'
out_dir='../processed_data/response_eQTL/residual/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

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
count=fread(count_fn)
count_df=as.data.frame(count[,2:ncol(count)])
rownames(count_df)=unlist(count[,1])


# Subset to protein coding and lncRNA:
count_df=subset_to_protein_and_lncRNA(count_df)


# Remove unexpressed genes: 
count_df=remove_unexpressed(count_df)


# Remove duplicate sample 021011: 
count_df$`021011.glucose`=NULL
count_df$`021011.galactose`=NULL


# Split count matrix into glucose and galactose: 
glu=count_df%>%select(contains('lucose'))
gal=count_df%>%select(contains('lactose'))


# Log transform:
glu_log=log2(glu+1)
gal_log=log2(gal+1)


# Remove covariates: 
glu_cov=fread(glu_cov_fn)
glu_res=matrix(NA,nrow=nrow(glu_log),ncol=ncol(glu_log))
rownames(glu_res)=rownames(glu_log)
colnames(glu_res)=colnames(glu_log)
for (i in 1:nrow(glu_log)){
	fit=lm(unlist(glu_log[i,])~V1+V2+V3+V4+V5+V6+V7+V8,data=glu_cov)
	glu_res[i,]=fit$residuals
}
glu_res=as.data.frame(glu_res)
glu_res=as.data.table(glu_res,keep.rownames=TRUE)
setnames(glu_res,'rn','gene_id')
fwrite(glu_res,sprintf('%s/glucose_residual.tsv',out_dir),sep='\t')



gal_cov=fread(gal_cov_fn)
gal_res=matrix(NA,nrow=nrow(gal_log),ncol=ncol(gal_log))
rownames(gal_res)=rownames(gal_log)
colnames(gal_res)=colnames(gal_log)
for (i in 1:nrow(gal_log)){
	fit=lm(unlist(gal_log[i,])~V1+V2+V3+V4+V5+V6+V7+V8+V9,data=gal_cov)
	gal_res[i,]=fit$residuals
}
gal_res=as.data.frame(gal_res)
gal_res=as.data.table(gal_res,keep.rownames=TRUE)
setnames(gal_res,'rn','gene_id')
fwrite(gal_res,sprintf('%s/galactose_residual.tsv',out_dir),sep='\t')


