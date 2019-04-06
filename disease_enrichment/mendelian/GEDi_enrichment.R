library(data.table)
library(cowplot)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(15)
source('utils/gtex_tissue_info.R')

tissue_fn = '../processed_data/rpe_specific_genes/rpe_specific_genes.GTExV7/tissue_kept.txt'
GEDi_fn = '../data/eye_disease/GEDi.txt'
epilepsy_fn = '../data/eye_disease/invitae_epilepsy.txt'
zscore_dir = '../processed_data/rpe_specific_genes/rpe_specific_genes.GTExV7/'
fig_dir = '../figures/disease_enrichment/mendelian/GEDi_enrichment/'
out_dir = '../processed_data/disease_enrichment/mendelian/GEDi_enrichment/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

read_tissue = function(tissue_fn){
	tissue = unname(unlist(fread(tissue_fn,header=FALSE, sep = '\t')))
	return(tissue)
}

read_GEDi = function(GEDi_fn){
	GEDi = fread(GEDi_fn,header=TRUE)
	return(GEDi)
}


read_zscore = function(zscore_fn){
	zscore = fread(zscore_fn)
	return(zscore)
}

calculate_t_stat = function(zscore_dir,tissue,disease_genes){
	t_stat = foreach(t = tissue,.combine = 'rbind')%dopar%{
		t_ = str_replace_all(t,' ','_')
		zscore_fn = sprintf('%s/all_genes.%s.txt',zscore_dir,t_)
		zscore = read_zscore(zscore_fn)
		zscore = zscore[type=='protein_coding'&chr%in%paste0('chr',1:22)]

		zscore[,GEDi := gene_name %in% disease_genes$gene_name]
		result = t.test(zscore~GEDi,data=zscore)
		x = with(result,data.frame(beta = estimate[2]-estimate[1],lb = -conf.int[2],ub = -conf.int[1],pval = p.value))
		x$tissue = t
		return(x)
	}
	rownames(t_stat) = NULL
	setDT(t_stat)
	setorder(t_stat,beta)
	t_stat[,tissue:=tissue_to_tissue_id[tissue]]
	t_stat[,tissue:=factor(tissue,levels=tissue)]
	return(t_stat)
}

plot_t_stat = function(t_stat){
	p = ggplot(t_stat,aes(x=tissue,y=beta,ymax=ub,ymin=lb,color=tissue))+
		geom_hline(yintercept = 0) + 
		geom_linerange(color='gray',size=1.5) +
		geom_point(size=3) + 
		xlab('') +
		ylab('Beta (95% CI)') + 
		scale_color_manual(values = tissue_color,guide='none') + 
		scale_x_discrete(labels = tissue_abbreviation) + 
		coord_flip()
	return(p)
}

tissue = read_tissue(tissue_fn)
tissue = c(tissue,'RPE')

GEDi = read_GEDi(GEDi_fn)
t_stat = calculate_t_stat(zscore_dir,tissue,GEDi)
fwrite(t_stat, '../processed_data/disease_enrichment/mendelian/GEDi_enrichment/t_stat.txt', sep = '\t')
t_stat[tissue=='RPE',pval] # 1.554947e-10
p = plot_t_stat(t_stat)
out_fn = sprintf('%s/GEDi_enrichment.rds',out_dir)
saveRDS(p,out_fn)
fig_fn = sprintf('%s/GEDi_enrichment.pdf',fig_dir)
save_plot(fig_fn,p,base_height=5.5,base_width=5.5)

disease_genes = read_GEDi(epilepsy_fn)
t_stat = calculate_t_stat(zscore_dir,tissue,disease_genes)
p = plot_t_stat(t_stat)
out_fn = sprintf('%s/epilepsy_enrichment.rds',out_dir)
saveRDS(p,out_fn)
fig_fn = sprintf('%s/epilepsy_enrichment.pdf',fig_dir)
save_plot(fig_fn,p,base_height=5.5,base_width=5.5)
