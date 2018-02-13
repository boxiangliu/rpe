#!/usr/bin/env Rscript
# boxiang liu
# durga
# perform multidimentional scaling 

# library:
library(dplyr)
library(data.table)
library('MASS')
library(cowplot)
library(stringr)
library(ggrepel)
library(foreach)
library(doMC)
registerDoMC(10)
library(preprocessCore)

# command line arguments: 
rpkm_file='../processed_data/mds/preprocess.all_tissue/combined.rpkm'
coldata_file='../processed_data/mds/preprocess.all_tissue/combined.col'
out_dir='../processed_data/rpe_specific_genes/'
if (!dir.exists(out_dir)){dir.create(out_dir)}

select_treatment=function(rpkm,col_data,treatment,format=c('data.table','matrix')){
	if ('Name'%in%colnames(rpkm)){
		gene_id=rpkm[,list(gene_id=Name)]
	}
	
	sample=col_data[tissue==treatment,sample]
	x=rpkm[,sample,with=FALSE]
	
	format=match.arg(format)
	if (format=='data.table'){
		cbind(gene_id,x)
	} else {
		x=as.matrix(x)
		rownames(x)=gene_id$gene_id
		return(x)
	}
}

subset_to_protein_and_lncRNA=function(x){
	gencode=fread('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed')[,5:6,with=F]
	setnames(gencode,c('gene_id','type'))
	y=x[rownames(x)%in%gencode[type=='lincRNA'|type=='protein_coding',gene_id],]
	return(y)
}

make_plot=function(cor,tissue=c('RPE (glu)','RPE (gal)'),color){
	x=data.table(cor=cor[colnames(cor)==tissue],tissue=colnames(cor))
	setorder(x,cor)
	x[,tissue:=factor(tissue,level=tissue)]

	p=ggplot(x,aes(cor,tissue,color=tissue))+geom_point()+scale_color_manual(values=color,guide='none')
	return(p)
}
 
main=function(){
	rpkm=fread(rpkm_file,header=T)

	rpe_rpkm=rpkm[,str_detect(colnames(rpkm),'(glu|gal)'),with=FALSE]
	glu_median=apply(rpkm[,str_detect(colnames(rpkm),'glucose'),with=FALSE],1,median)
	gal_median=apply(rpkm[,str_detect(colnames(rpkm),'galactose'),with=FALSE],1,median)
	cor(glu_median,gal_median) # 0.9964466

	col_data=fread(coldata_file,header=T)
	col_data[,tissue:=str_replace(tissue,'gal','')]
	col_data[,tissue:=str_replace(tissue,'glu','')]

	median=foreach(tissue=unique(col_data$tissue),.combine='cbind')%dopar%{
		x=select_treatment(rpkm,col_data,tissue,'matrix')
		x_median=apply(x,1,median)
		y=data.table(x_median)
		setnames(y,'x_median',tissue)
		return(y)
	}
	median=as.matrix(median)
	rownames(median)=rpkm$Name

	median_filt=subset_to_protein_and_lncRNA(median)

	corr=cor(median_filt)
	diag(corr)=0 # set diagonal to 0 to keep current tissue in each iteration.

	median_filt=median_filt[,c(ncol(median_filt),1:(ncol(median_filt)-1))]

	threshold=0.96
	n_tissue_kept=0
	n_tissue_remaining=nrow(corr)
	neighboring_tissue=list()
	while(n_tissue_remaining>0){
		n_tissue_kept=n_tissue_kept+1

		tissue_to_remove=names(which(corr[n_tissue_kept,]>=threshold))
		neighboring_tissue[[rownames(corr)[n_tissue_kept]]]=tissue_to_remove

		tissue_to_keep=which(corr[n_tissue_kept,]<threshold)
		corr=corr[tissue_to_keep,tissue_to_keep]
		
		n_tissue_remaining=nrow(corr)-n_tissue_kept
	}

	tissue_kept=colnames(corr)
	write.table(tissue_kept[1:(length(tissue_kept)-1)],paste0(out_dir,'/tissue_kept.txt'),sep='\t',row.names=FALSE,col.names=FALSE)

	median_filt_independent=median_filt[,tissue_kept]

	sd=apply(median_filt_independent,1,sd)
	mean=apply(median_filt_independent,1,mean)
	idx=which((median_filt_independent[,'RPE ()']-mean)/sd>4)

	out=data.table(gene_id=rownames(median_filt_independent),
		rpkm=median_filt_independent[,'RPE ()'],
		mean=mean,
		sd=sd)
	out[,zscore:=(rpkm-mean)/sd]
	gencode=fread('../data/gtex/gencode.v19.genes.v6p.hg19.bed',col.names=c('chr','start','stop','strand','gene_id','gene_name','type'))
	out=merge(out,gencode,by='gene_id')

	fwrite(out,paste0(out_dir,'/all_genes.txt'),sep='\t')
	fwrite(out[zscore>4,],paste0(out_dir,'/rpe_specific_genes.txt'),sep='\t')
}

main()
