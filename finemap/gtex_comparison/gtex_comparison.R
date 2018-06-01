library(data.table)
library(stringr)
library(cowplot)
library(ggrepel)
source('utils/gtex_tissue_info.R')
source('utils/genome_annotation.R')
library(foreach)

myopia_gtex_fn='../processed_data/finemap/manhattan/2018-01-31_11-50-44_rpe_23andme_gtex/23andme_myopia_prepared_txt_gz_finemap_clpp_status.txt'
myopia_rpe_fn='../processed_data/finemap/manhattan/2018-02-27_11-56-08_rpe_23andme_updated/23andme_myopia_prepared_txt_gz_finemap_clpp_status.txt'
amd_rpe_fn = '../processed_data/finemap/manhattan/2018-02-27_19-03-38_rpe_amd_updated/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'
amd_gtex_fn = '../processed_data/finemap/manhattan/2018-02-08_16-09-37_rpe_amd_gtex/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'
fig_dir='../figures/finemap/gtex_comparison/gtex_comparison/'

if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

read_finemap=function(fm_fn,threshold=5e-5,subset=NULL){
	tmp=fread(fm_fn)
	if (ncol(tmp)==8){
		setnames(tmp,c('snp','eqtl_file','gwas_file','gene_id','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
	} else if (ncol(tmp)==9){
		setnames(tmp,c('snp','eqtl_file','gwas_file','gene_id','conditional_level','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
	} else {
		stop('Check input column names!')
	}
	tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
	tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
	fm=tmp[,list(chrom,pos,tissue=eqtl_file,gene_id,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR')]

	set.seed(42)
	fm[,rank:=rank(-y,ties.method='random'),by=c('gene_id','tissue')]
	fm=fm[rank==1]
	fm$rank=NULL
	fm=fm[gwas_logp> -log10(threshold)|eqtl_logp> -log10(threshold)]
	fm[,tissue:=gsub('_Analysis_cis_eqtl_gz|_eqtls_txt_gz|_allpairs_txt_gz','',tissue)]
	if (!is.null(subset)) {
		fm=fm[gene_id%in%subset]
	}
	return(fm)
}

munge_rpe = function(rpe){
	x = copy(rpe)
	x = x[str_detect(tissue,'(^glucose$|^galactose$)')]
	x = x[,tissue := paste0('RPE - ',tissue)]
	x[,gene_id := str_split_fixed(gene_id,'_',2)[,1]]

	return(x)
}

plot_clpp=function(data,color_map,top=length(unique(data$tissue)),label_top=length(unique(data$tissue))){
	setorder(data,-y)
	data=data[1:top]
	setorder(data,y)
	data[,tissue:=factor(tissue,tissue)]
	p=ggplot(data,aes(y,tissue,color=tissue))+
		geom_point()+
		scale_color_manual(values=tissue_color,guide='none')+
		scale_y_discrete(labels = tissue_abbreviation) + 
		xlab(sprintf('%s CLPP',unique(data$gene_name)))+
		ylab('')+ 
		theme(axis.text.y=element_text(color=ifelse(str_detect(data$tissue,'RPE'),'purple','black')))
	return(p)
}



# AMD:
amd_gtex = read_finemap(amd_gtex_fn,threshold=1)
amd_rpe = read_finemap(amd_rpe_fn,threshold=1)
amd_rpe = munge_rpe(amd_rpe)
amd_gtex = merge(amd_gtex,gene_annotation[,list(gene_id,gene_name)],by='gene_id')
amd_rpe = merge(amd_rpe,gene_annotation[,list(gene_id,gene_name)],by='gene_id')
amd=rbind(amd_gtex,amd_rpe)

gene_list = c('RDH5','PARP12')
p_list = foreach(i = gene_list)%do%{
	data=amd[gene_name==i]
	p=plot_clpp(data,color_map,top=30)
	return(p)
}
p = plot_grid(plotlist=p_list,nrow=1,labels=LETTERS[1:2])
save_plot(sprintf('%s/gtex_comparison_amd.pdf',fig_dir),p,base_height=5,base_width=8.5)


# Myopia:
myopia_gtex = read_finemap(myopia_gtex_fn,threshold=1)
myopia_rpe = read_finemap(myopia_rpe_fn,threshold=1)
myopia_rpe = munge_rpe(myopia_rpe)
gene_annotation = read_gencode(gencode_fn)
myopia_gtex = merge(myopia_gtex,gene_annotation[,list(gene_id,gene_name)],by='gene_id')
myopia_rpe = merge(myopia_rpe,gene_annotation[,list(gene_id,gene_name)],by='gene_id')

myopia=rbind(myopia_gtex,myopia_rpe)
gene_list = c('RDH5','CLU','APH1B','BCKDHA')
p_list = foreach(i = gene_list)%do%{
	data=myopia[gene_name==i]
	p=plot_clpp(data,color_map,top=30)
	return(p)
}
p = plot_grid(plotlist=p_list,nrow=2,labels=LETTERS[1:4])
save_plot(sprintf('%s/gtex_comparison_myopia.pdf',fig_dir),p,base_height=10,base_width=8)

