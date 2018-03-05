library(data.table)
library(cowplot)
library(ggrepel)
library(manhattan)
library(locuscomparer)
library(stringr)

clpp_fn='../processed_data/finemap/manhattan/2018-02-02_12-36-28_rpe_amd/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'
fig_dir='../figures/finemap/manhattan/manhattan_amd/'
out_dir='../processed_data/finemap/manhattan/manhattan_amd/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_finemap=function(fm_fn,threshold=5e-5){
	if (grepl('clpp_status',fm_fn)) {
		tmp=fread(fm_fn,col.names=c('snp','eqtl_file','gwas_file','gene_name','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
		tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
		tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
		tmp[,condition:='']
		tmp[,condition:=ifelse(str_detect(eqtl_file,'glucose'),'Glucose',condition)]
		tmp[,condition:=ifelse(str_detect(eqtl_file,'galactose'),'Galactose',condition)]
		fm=tmp[grepl('eqtl',eqtl_file),list(chrom,pos,gene_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR',condition)]

	} else {
		fm=fread(fm_fn)[,list(chrom=chrom,pos,gene_name=gene,y=clpp_score,gwas_logp=gwas_log_pval,eqtl_logp=eqtl_log_pval,method='eCAVIAR')]
	}

	set.seed(42)
	fm=fm[gwas_logp> -log10(threshold)|eqtl_logp> -log10(threshold)]
	fm[,gene_id:=str_split_fixed(gene_name,'_',2)[,1]]
	fm[,gene_name:=str_split_fixed(gene_name,'_',2)[,2]]
	return(fm)
}


data=read_finemap(clpp_fn,threshold=1)
cutoff = 0.05
data[,condition:=factor(condition,level=c('Glucose','Galactose'))]
data[,label:=ifelse(y>cutoff,gene_name,'')]
data[,chrom:=paste0('chr',chrom)]
data[condition=='Glucose',y:=-y]
dummy=data.table(condition=c('Glucose','Galactose'),y=c(-cutoff,cutoff))

p=manhattan(data,build='hg19')+
	facet_grid(condition~.,scale='free_y')+
	scale_y_continuous(labels=function(x){abs(x)})+
	geom_text_repel(aes(label=label))+
	geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
	ylab(paste('CLPP'))


saveRDS(list(data,dummy,p),sprintf('%s/manhattan.rds',out_dir))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=8,base_height=4)