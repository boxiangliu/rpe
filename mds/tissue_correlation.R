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

# command line arguments: 
rpkm_file='../processed_data/mds/preprocess.all_tissue/combined.rpkm'
coldata_file='../processed_data/mds/preprocess.all_tissue/combined.col'
fig_dir='../figures/mds/tissue_correlation/'
if (!dir.exists(fig_dir)){dir.create(fig_dir)}

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

# read input: 
main=function(){
	rpkm=fread(rpkm_file,header=T)
	col_data=fread(coldata_file,header=T)
	median=foreach(tissue=unique(col_data$tissue),.combine='cbind')%dopar%{
		x=select_treatment(rpkm,col_data,tissue,'matrix')
		x_median=apply(x,1,median)
		y=data.table(x_median)
		setnames(y,'x_median',tissue)
		return(y)
	}
	rownames(median)=rpkm$Name

	median_filt=subset_to_protein_and_lncRNA(median)

	cor=cor(median_filt)
	gtex_color=fread('/srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gtex_tissue_colors.with_hcasmc.txt')[,list(tissue_site_detail,tissue_color_hex)]
	gtex_color=rbind(gtex_color,data.table(tissue_site_detail=c('RPE (glu)','RPE (gal)'),tissue_color_hex=c('#000000','#000000')))
	color=gtex_color$tissue_color_hex
	names(color)=gtex_color$tissue_site_detail

	pdf(sprintf('%s/tissue_correlation_pdf',fig_dir))
	print(make_plot(cor,'RPE (gal)',color))
	print(make_plot(cor,'RPE (glu)',color))
	dev.off()



}
