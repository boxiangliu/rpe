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
library(manhattan)
library(ggrepel)

# command line arguments: 
rpkm_file='../processed_data/mds/preprocess.GTExV7/combined_logxp2.rpkm'
coldata_file='../processed_data/mds/preprocess.GTExV7/combined.col'
out_dir='../processed_data/rpe_specific_genes/rpe_specific_genes.GTExV7/'
fig_dir='../figures/rpe_specific_genes/rpe_specific_genes.GTExV7/'
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
 
calculate_median_rpkm = function(rpkm,col_data){
	median=foreach(tissue=unique(col_data$tissue),.combine='cbind')%dopar%{
		x=select_treatment(rpkm,col_data,tissue,'matrix')
		x_median=apply(x,1,median)
		y=data.table(x_median)
		setnames(y,'x_median',tissue)
		return(y)
	}
	median=as.matrix(median)
	rownames(median)=rpkm$Name
	return(median)
}

select_indepedent_tissue = function(corr,out_dir,threshold=0.96){
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
	return(tissue_kept)
}

make_plot_data = function(out,biotype){
	plot_data = copy(out)
	plot_data[,c('chrom','pos','y'):=list(chr,start,zscore)]
	plot_data = plot_data[!chr%in%c('chrX','chrY','chrM'),]
	plot_data[,color:=ifelse(abs(zscore)>4,'red',NA)]
	plot_data[color=='red',color:=ifelse(chr%in%paste0('chr',seq(1,22,2)),'red','pink')]
	mhc_genes = plot_data[chr=='chr6'&start>28477797&stop<33448354,gene_id]
	plot_data = plot_data[!gene_id%in%mhc_genes,] # remove MHC genes
	plot_data = plot_data[type==biotype]
	plot_data[,rank := rank(-zscore)]
	plot_data[,label := ifelse(rank<=15,gene_name,NA)]
	plot_data = add_cumulative_pos(plot_data, 'hg19')
	return(plot_data)
}

plot_zscore = function(plot_data,top){
	p = manhattan(plot_data,build = 'hg19') + 
		geom_text_repel(
			data = plot_data[rank<=top],
			aes(label = gene_name),
			color='black',
			direction='x',
			segment.color = "grey",
			nudge_y=9-plot_data[rank<=top]$y,
			angle=90)+
		ylab('RPE specificity (z-score)')+
		scale_y_continuous(breaks = seq(-2,9,2),limits = c(-3,9))
	return(p)
}

# Read RPKM:
rpkm=fread(rpkm_file,header=T)

# Calculate glucose and galatose correlation:
rpe_rpkm=rpkm[,str_detect(colnames(rpkm),'(lucose|alactose)'),with=FALSE]
glu_median=apply(rpkm[,str_detect(colnames(rpkm),'lucose'),with=FALSE],1,median)
gal_median=apply(rpkm[,str_detect(colnames(rpkm),'alactose'),with=FALSE],1,median)
cor(glu_median,gal_median) # 0.9975472

# Read column data:
col_data=fread(coldata_file,header=T)
col_data[,tissue:=trimws(str_replace(tissue,'\\(gal\\)',''))]
col_data[,tissue:=trimws(str_replace(tissue,'\\(glu\\)',''))]

# Calculate median RPKM:
median = calculate_median_rpkm(rpkm,col_data)
median_filt=subset_to_protein_and_lncRNA(median)

# Select independent tissues:
corr=cor(median_filt)
diag(corr)=0 # set diagonal to 0 to keep current tissue in each iteration.

median_filt=median_filt[,c(ncol(median_filt),1:(ncol(median_filt)-1))] # move RPE to first column


tissue_kept = select_indepedent_tissue(corr,out_dir)
median_filt_independent=median_filt[,tissue_kept]

# Calculate mean and sd:
sd=apply(median_filt_independent,1,sd)
mean=apply(median_filt_independent,1,mean)

# Calculate z-score for every tissue:
foreach(tissue = colnames(median_filt_independent))%dopar%{
	print(tissue)
	out=data.table(gene_id=rownames(median_filt_independent),
		rpkm=median_filt_independent[,tissue],
		mean=mean,
		sd=sd)
	out[,zscore:=(rpkm-mean)/sd]
	gencode=fread('../data/gtex/gencode.v19.genes.v6p.hg19.bed',col.names=c('chr','start','stop','strand','gene_id','gene_name','type'))
	out=merge(out,gencode,by='gene_id')

	tissue = str_replace_all(tissue,' ','_')
	fwrite(out,sprintf('%s/all_genes.%s.txt',out_dir,tissue),sep='\t')
	fwrite(out[zscore>4,],sprintf('%s/specific_genes.%s.txt',out_dir,tissue),sep='\t')
}


# Caculate again for RPE for plotting:
out=data.table(gene_id=rownames(median_filt_independent),
	rpkm=median_filt_independent[,'RPE'],
	mean=mean,
	sd=sd)
out[,zscore:=(rpkm-mean)/sd]
gencode=fread('../data/gtex/gencode.v19.genes.v6p.hg19.bed',col.names=c('chr','start','stop','strand','gene_id','gene_name','type'))
out=merge(out,gencode,by='gene_id')


# Plot z-score manhattan plot:
plot_data = make_plot_data(out,'protein_coding')
p = plot_zscore(plot_data,top=17)
save_plot(sprintf('%s/manhattan_zscore_proteinCoding.pdf',fig_dir),p,base_width=8,base_height=4)

plot_data = make_plot_data(out,'lincRNA')
p = plot_zscore(plot_data,top=10)
save_plot(sprintf('%s/manhattan_zscore_lincRNA.pdf',fig_dir),p,base_width=8,base_height=4.5)