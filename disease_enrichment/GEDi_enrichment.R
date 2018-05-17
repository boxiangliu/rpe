library(data.table)
library(cowplot)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(15)

tissue_fn = '../processed_data/rpe_specific_genes.GTExV7/tissue_kept.txt'
GEDi_fn = '../data/eye_disease/GEDi.txt'
zscore_fn = '../processed_data/rpe_specific_genes.GTExV7/all_genes.RPE.txt'
zscore_dir = '../processed_data/rpe_specific_genes.GTExV7/'
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'

read_tissue = function(tissue_fn){
	tissue = unname(unlist(fread(tissue_fn,header=FALSE)))
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

read_tissue_color = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_color = gtex_tissue_color$tissue_color_hex
	names(tissue_color) = gtex_tissue_color$tissue_site_detail_id
	tissue_color = c(tissue_color,c(`RPE - glucose` = '#FF0000', `RPE - galactose` = '#00FF00',RPE = '#FF0000'))
	return(tissue_color)
}

read_tissue_abbreviation = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_abbreviation = gtex_tissue_color$tissue_abbreviation
	names(tissue_abbreviation) = gtex_tissue_color$tissue_site_detail_id
	tissue_abbreviation = c(tissue_abbreviation,c(`RPE - glucose` = 'RPE - glucose', `RPE - galactose` = 'RPE - galactose'))
	return(tissue_abbreviation)
}

get_tissue_to_tissue_id = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	tissue_to_tissue_id = gtex_tissue_color$tissue_site_detail_id
	names(tissue_to_tissue_id) = gtex_tissue_color$tissue_site_detail
	tissue_to_tissue_id = c(tissue_to_tissue_id,`RPE (glu)`="RPE - glucose",`RPE (gal)`="RPE - galactose",RPE='RPE')
	return(tissue_to_tissue_id)
}

tissue = read_tissue(tissue_fn)
tissue = c(tissue,'RPE')
GEDi = read_GEDi(GEDi_fn)

t_stat = foreach(t = tissue,.combine = 'rbind')%dopar%{
	t_ = str_replace_all(t,' ','_')
	zscore_fn = sprintf('%s/all_genes.%s.txt',zscore_dir,t_)
	zscore = read_zscore(zscore_fn)
	zscore = zscore[type=='protein_coding'&chr%in%paste0('chr',1:22)]

	zscore[,GEDi := gene_name %in% GEDi$gene_name]
	result = t.test(zscore~GEDi,data=zscore)
	x = with(result,data.frame(beta = estimate[2]-estimate[1],lb = -conf.int[2],ub = -conf.int[1],pval = p.value))
	x$tissue = t
	return(x)
}

tissue_to_tissue_id = get_tissue_to_tissue_id(gtex_tissue_color_fn)
tissue_color = read_tissue_color(gtex_tissue_color_fn)
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)

rownames(t_stat) = NULL
setDT(t_stat)
setorder(t_stat,beta)
t_stat[,tissue:=tissue_to_tissue_id[tissue]]
t_stat[,tissue:=factor(tissue,levels=tissue)]

ggplot(t_stat,aes(x=tissue,y=beta,ymax=ub,ymin=lb,color=tissue))+
geom_hline(yintercept = 0) + 
geom_linerange(color='gray',size=1.5) +
geom_point(size=3) + 
xlab('') +
ylab('Enrichment beta (95% CI)') + 
scale_color_manual(values = tissue_color,guide='none') + 
scale_x_discrete(labels = tissue_abbreviation) + 
coord_flip() 

ggplot(zscore,aes(GEDi,zscore)) + geom_violin()
ggplot(zscore,aes(GEDi,zscore)) + geom_boxplot()
ggplot(zscore,aes(zscore,color=GEDi)) + geom_density()
