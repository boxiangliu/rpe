library(data.table)
library(cowplot)
library(ggrepel)
library(locuscomparer)
library(leafcutter)
library(stringr)
source('finemap/manhattan/utils.R')

clpp_fn='../processed_data/finemap/manhattan/2018-06-06_15-06-09_rpe_23andme/23andme_myopia_prepared_txt_gz_finemap_clpp_status.txt'
exon_file='/srv/persistent/bliu2/tools/leafcutter/leafcutter/data/gencode19_exons.txt.gz'
fig_dir='../figures/finemap/manhattan/manhattan_myopia_sqtl/'
out_dir='../processed_data/finemap/manhattan/manhattan_myopia_sqtl/'
if(!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}
if(!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

read_finemap=function(fm_fn, threshold1, threshold2){
	if (grepl('clpp_status',fm_fn)) {
		tmp=fread(fm_fn,select=1:8,col.names=c('snp','eqtl_file','gwas_file','clu_name','n_tested_snps','clpp_score','gwas_logp','eqtl_logp'))
		tmp[,chrom:=str_split_fixed(snp,'_',2)[,1]]
		tmp[,pos:=as.integer(str_split_fixed(snp,'_',2)[,2])]
		tmp[,condition:='']
		tmp[,condition:=ifelse(str_detect(eqtl_file,'glucose'),'Glucose',condition)]
		tmp[,condition:=ifelse(str_detect(eqtl_file,'galactose'),'Galactose',condition)]
		fm=tmp[grepl('sqtl',eqtl_file),list(chrom,pos,clu_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR',condition)]

	} else {
		fm=fread(fm_fn)[,list(chrom=chrom,pos,clu_name=gene,y=clpp_score,gwas_logp=gwas_log_pval,eqtl_logp=eqtl_log_pval,method='eCAVIAR')]
	}

	set.seed(42)
	fm=fm[ eqtl_logp> -log10(threshold1) & gwas_logp> -log10(threshold2)]
	return(fm)
}

assign_gene=function(exon_file,clusters){
	exons_table     = read.table(exon_file, header=T, stringsAsFactors = F)
	intron_meta     = get_intron_meta(clusters)
	exons_table$chr = add_chr(exons_table$chr)
	intron_meta$chr = add_chr(intron_meta$chr)
	clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
	clu_map = data.frame(clusters,clu=gsub(':[0-9]+:[0-9]+','',clusters))
	clu_gene_map = merge(clu_gene_map,clu_map,by='clu')
	map = clu_gene_map$genes
	names(map) = clu_gene_map$clusters
	return(map[clusters])
}

data=read_finemap(clpp_fn,threshold1=1e-5,threshold2=1e-4)
data$clu_name=paste0('chr',gsub('clu:','clu_',gsub('\\.',':',data$clu_name)))
data$gene_name=assign_gene(exon_file,data$clu_name)
data[,rank:=rank(-y,ties.method='first'),by=c('gene_name','condition')]
data = data[rank==1]
data$rank = NULL
cutoff = 0.01
data[,condition:=factor(condition,level=c('Glucose','Galactose'))]
data[,label:=ifelse(y>cutoff,gene_name,'')]
data[,chrom:=paste0('chr',chrom)]
p = plot_manhattan(data,cutoff=cutoff)

saveRDS(list(data,p),sprintf('%s/manhattan.rds',out_dir))
save_plot(sprintf('%s/manhattan.pdf',fig_dir),p,base_width=6,base_height=4)