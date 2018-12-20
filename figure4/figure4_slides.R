# library(cowplot)
library(ggplot2)
library(data.table)
library(foreach)
source('utils/gtex_tissue_info.R')

response_eQTL_dir = '../processed_data/QTL_landscape/response_eQTL/'
glucose_eQTL_rds = paste0(response_eQTL_dir,'glucose_specific_eQTL_plots.rds')
galactose_eQTL_rds = paste0(response_eQTL_dir,'galactose_specific_eQTL_plots.rds')
shared_eQTL_rds = paste0(response_eQTL_dir,'shared_eQTL_plots.rds')
rpe_specific_eQTL_rds = '../processed_data/rpe_specific_eQTL/specific_eGenes_v2/rpe_specific_eGenes.rds'
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
fig_dir = '../figures/figure4/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


plot_treatment_specific_eQTL = function(rds, genes, labels){
	eQTL = readRDS(rds)
	p = foreach(i = 1:length(genes))%do%{
		p = eQTL[[genes[i]]] + 
		scale_fill_manual(values=c(Glucose='#F87660',Galactose='#619CFF'),guide='none') + 
		theme(strip.text = element_blank(),strip.background=element_blank()) + 
		ylab('')
		return(p)
	}
	plot_grid(plotlist = p,nrow = 1, align = 'h', labels = labels) 
}

glucose_eQTL_plot = plot_treatment_specific_eQTL(glucose_eQTL_rds,genes = c('ABCA1'),labels = '')