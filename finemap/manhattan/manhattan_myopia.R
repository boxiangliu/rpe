library(data.table)
library(cowplot)
library(ggrepel)
library(locuscomparer)
library(stringr)
source('finemap/manhattan/utils.R')

glu_fn='../processed_data/finemap/manhattan/2018-01-26_12-23-26_rpe_23andme/manhattan/23andme_myopia_prepared_txt_gz_glucose_eqtls_txt_gz.txt'
gal_fn='../processed_data/finemap/manhattan/2018-01-26_12-23-26_rpe_23andme/manhattan/23andme_myopia_prepared_txt_gz_galactose_eqtls_txt_gz.txt'
fig_dir='../figures/finemap/manhattan/manhattan_myopia/'
out_dir='../processed_data/finemap/manhattan/manhattan_myopia/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_finemap=function(fm_fn,threshold=5e-5){
	if (grepl('clpp_status',fm_fn)) {
		tmp=fread(fm_fn,col.names=c('snp','eqtl_file','gwas_file','gene_name','conditional_level','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
		tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
		tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
		fm=tmp[grepl('eqtl',eqtl_file),list(chrom,pos,gene_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR')]
	} else {
		fm=fread(fm_fn)[,list(chrom=chrom,pos,gene_name=gene,y=clpp_score,gwas_logp=gwas_log_pval,eqtl_logp=eqtl_log_pval,method='eCAVIAR')]
	}

	set.seed(42)
	fm[,rank:=rank(-y,ties.method='random'),by='gene_name']
	fm=fm[rank==1]
	fm$rank=NULL
	fm=fm[gwas_logp> -log10(threshold)|eqtl_logp> -log10(threshold)]
	fm[,gene_id:=str_split_fixed(gene_name,'_',2)[,1]]
	fm[,gene_name:=str_split_fixed(gene_name,'_',2)[,2]]
	return(fm)
}


glu=read_finemap(glu_fn)[,condition:='Glucose']
gal=read_finemap(gal_fn)[,condition:='Galactose']

cutoff = 0.05
data=rbind(glu,gal)
data[,condition:=factor(condition,level=c('Glucose','Galactose'))]
data[,label:=ifelse(y>cutoff,gene_name,'')]
data[,chrom:=paste0('chr',chrom)]
p = plot_manhattan(data,cutoff=cutoff)

saveRDS(list(data,p),sprintf('%s/manhattan.rds',out_dir))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=8,base_height=4)


