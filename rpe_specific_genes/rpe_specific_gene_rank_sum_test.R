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
out_dir='../processed_data/rpe_specific_genes/rpe_specific_genes_rank_sum_test/'
fig_dir='../figures/rpe_specific_genes/rpe_specific_genes_rank_sum_test/'
if (!dir.exists(out_dir)){dir.create(out_dir)}
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
 
main=function(){
	rpkm=fread(rpkm_file,header=T)
	col_data=fread(coldata_file,header=T)
	tissue_kept=unlist(fread('../processed_data/rpe_specific_genes/tissue_kept.txt',header=FALSE))

	gtex_col_data=col_data[!str_detect(tissue,'RPE')&tissue%in%tissue_kept,]

	# Select non-overlapping individuals for each tissue:
	gtex_col_data[,individual:=str_extract(sample,'GTEX-[A-Z0-9]+?(?=-)')]
	sample_size=gtex_col_data[,list(size=.N),by=tissue]
	setorder(sample_size,size)
	selected_indv=list()
	selected_gtex_col_data=data.table()
	for (tiss in sample_size$tissue){
		print(tiss)
		indv=gtex_col_data[tissue==tiss&!(individual%in%selected_indv),individual]
		if (length(indv)>24){
			indv=sample(indv,24)
		}
		selected_indv=c(selected_indv,indv)
		selected_gtex_col_data=rbind(selected_gtex_col_data,gtex_col_data[tissue==tiss&individual%in%indv])
	}
	fwrite(selected_gtex_col_data[,list(size=.N),by='tissue'],paste0(out_dir,'tissue_sample_size.txt'),sep='\t')

	# Rank-sum test:
	gtex_rpkm=rpkm[,selected_gtex_col_data$sample,with=FALSE]
	glucose_rpkm=rpkm[,str_detect(colnames(rpkm),'glucose'),with=FALSE]
	combined_rpkm=as.data.frame(cbind(gtex_rpkm,glucose_rpkm))
	rownames(combined_rpkm)=rpkm$Name

	qnorm_rpkm=normalize.quantiles(as.matrix(combined_rpkm),copy=TRUE)
	rownames(qnorm_rpkm)=rownames(combined_rpkm)
	colnames(qnorm_rpkm)=colnames(combined_rpkm)
	gtex_rpkm=qnorm_rpkm[,selected_gtex_col_data$sample]
	glucose_rpkm=qnorm_rpkm[,str_detect(colnames(rpkm),'glucose')]
	pval=foreach(i=1:nrow(gtex_rpkm),.combine='c')%dopar%{
		wilcox.test(unlist(gtex_rpkm[i,]),unlist(glucose_rpkm[i,]),alternative='less')$p.value
	}

	names(pval)=rpkm$Name
	pval=data.frame(gene_id=rpkm$Name,pval)
	pval_protein_and_lncRNA=subset_to_protein_and_lncRNA(pval)


	pdf(paste0(fig_dir,'pval_histogram.pdf'))
	hist(pval_protein_and_lncRNA$pval,breaks=1000)
	dev.off()

	sum(pval_protein_and_lncRNA$pval<0.05/length(pval_protein_and_lncRNA)) # 1776

}


main()
