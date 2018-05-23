library(data.table)
library(cowplot)
library(stringr)

source('utils/gtex_tissue_info.R')

amd_fn='../processed_data/disease_enrichment/ld_score_regression/partition_heritability_merged/top500/amd.merged.nobaseline.results'
myopia_fn = '../processed_data/disease_enrichment/ld_score_regression/partition_heritability_merged/top500/myopia.merged.nobaseline.results'
fig_dir = '../figures/disease_enrichment/ld_score_regression/plot_heritability_enrichment/'
out_dir = '../processed_data/disease_enrichment/ld_score_regression/plot_heritability_enrichment/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_ldscore_regression = function(fn){
	x = fread(fn)
	x[,Coefficient_pval := 2*pnorm(abs(`Coefficient_z-score`),lower.tail=FALSE)]
	x[,Category:=str_replace(Category,'L2_0','')]
	x[,Category:=str_replace_all(Category,'_',' ')]
	x[,Category:=tissue_to_tissue_id[Category]]
	return(x)
}

plot_coefficient_pval = function(x){
	setorder(x,-Coefficient_pval)
	x[,Category := factor(Category,Category)]
	p = ggplot(x,aes(Category,-log10(Coefficient_pval),fill=Category)) + 
		geom_bar(stat='identity',color='black') + 
		scale_fill_manual(values = tissue_color,guide='none') + 
		xlab('') + 
		ylab(expression(-log[10]*'(coefficient p-value)')) + 
		scale_x_discrete(labels = tissue_abbreviation) + 
		coord_flip()
	return(p)
}

plot_enrichment_pval = function(x){
	setorder(x,-Enrichment_p)
	x[,Category := factor(Category,Category)]
	p = ggplot(x,aes(Category,-log10(Enrichment_p),fill=Category)) + 
		geom_bar(stat='identity',color='black') + 
		scale_fill_manual(values = tissue_color,guide='none') + 
		xlab('') + 
		ylab(expression(-log[10]*'(enrichment p-value)')) + 
		scale_x_discrete(labels = tissue_abbreviation) + 
		coord_flip()
	return(p)
}

amd = read_ldscore_regression(amd_fn)
myopia = read_ldscore_regression(myopia_fn)

p1=plot_coefficient_pval(amd)
fig_fn = sprintf('%s/amd_coefficient_pval.pdf',fig_dir)
save_plot(fig_fn,p1,base_height=6,base_width=6)
out_fn = sprintf('%s/amd_coefficient_pval.rds',out_dir)
saveRDS(p1,out_fn)

p2=plot_coefficient_pval(myopia)
fig_fn = sprintf('%s/myopia_coefficient_pval.pdf',fig_dir)
save_plot(fig_fn,p2,base_height=6,base_width=6)
out_fn = sprintf('%s/myopia_coefficient_pval.rds',out_dir)
saveRDS(p2,out_fn)

p3 = plot_enrichment_pval(amd)
fig_fn = sprintf('%s/amd_enrichment_pval.pdf',fig_dir)
save_plot(fig_fn,p3,base_height=6,base_width=6)

p4 = plot_enrichment_pval(myopia)
fig_fn = sprintf('%s/myopia_enrichment_pval.pdf',fig_dir)
save_plot(fig_fn,p4,base_height=6,base_width=6)

# Make supplementary figure:
p3 = p3 + ggtitle('AMD')
p4 = p4 + ggtitle('Myopia')
p = plot_grid(p3,p4,labels=c('A','B'))
fig_fn = sprintf('%s/grid_enrichment_pval.pdf',fig_dir)
save_plot(fig_fn,p,base_height=6,base_width=9)