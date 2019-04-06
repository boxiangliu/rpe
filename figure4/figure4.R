library(cowplot)
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
out_dir = '../processed_data/figure4/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

plot_blank = function(){
	ggplot() + geom_blank()
}

plot_treatment_specific_eQTL = function(rds, genes, labels){
	eQTL = readRDS(rds)
	p = foreach(i = 1:length(genes))%do%{
		p = eQTL[[genes[i]]] + 
		scale_fill_manual(values=c(Glucose='#F87660',Galactose='#619CFF'),guide='none') + 
		theme(strip.text = element_blank(),strip.background=element_blank()) + 
		ylab('Expr.')
		return(p)
	}
	plot_grid(plotlist = p,nrow = 1, align = 'h', labels = labels) 
}

make_plot_data = function(gene, GTEx,glucose_eGenes,galactose_eGenes){
	data1 = GTEx[gene_id == gene,list(tissue,FDR = qval)]
	data2 = glucose_eGenes[family == gene, list(tissue = "RPE - glucose",FDR = fam_p)]
	data3 = galactose_eGenes[family == gene, list(tissue = 'RPE - galactose', FDR = fam_p)]
	data = rbind(data1, data2, data3)
	return(data)
}

plot_FDR = function(data){
	data[,logFDR := -log10(FDR)]
	setorder(data,logFDR)
	data[,tissue := factor(tissue,tissue)]
	ggplot(data, aes(tissue, logFDR, fill = tissue)) + 
		geom_bar(stat = 'identity',color='black') + 
		geom_hline(yintercept = -log10(0.1),color='black',linetype='dashed') +
		scale_fill_manual(guide='none',values = tissue_color) + 
		scale_x_discrete(breaks = names(tissue_abbreviation), labels = tissue_abbreviation) + 
		ylab(expression('-'*log[10]*'(FDR)')) + 
		xlab('') + 
		coord_flip()
}

plot_rpe_specific_eQTL = function(rpe_specific_eQTL_rds,labels,out_dir = NULL){
	tissue_color = read_tissue_color(gtex_tissue_color_fn)
	tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)

	result = readRDS(rpe_specific_eQTL_rds)
	GTEx_eQTL = result[['GTEx_eQTL']]
	GTEx = result[['GTEx']]
	glucose_eGenes = result[['glucose_eGenes']]
	galactose_eGenes = result[['galactose_eGenes']]
	tissue_color = result[['tissue_color']]
	tissue_abbreviation = result[['tissue_abbreviation']]

	gene_list = unique(GTEx_eQTL[eQTL == FALSE,gene_id])
	p = foreach(gene = gene_list)%do%{
		data = make_plot_data(gene, GTEx, glucose_eGenes, galactose_eGenes)
		gene_name = GTEx_eQTL[gene_id==gene,unique(gene_name)]
		if (!is.null(out_dir)) {
			fwrite(data, sprintf('%s/%s.txt',out_dir,gene_name),sep='\t')
		}
		plot_FDR(data) + ggtitle(gene_name) + theme(axis.title=element_text(size=12))
	}

	col2 = plot_grid(p[[2]]+xlab(''),p[[3]],align='v',nrow=2,labels = labels[2:3],rel_heights = c(6,8))
	grid_p = plot_grid(p[[1]],col2,nrow=1,labels=c(labels[1],''))
	return(grid_p)
}


plot_entire_figure = function(labels = letters[1:7]){
	blank = plot_blank()
	glucose_eQTL_plot = plot_treatment_specific_eQTL(glucose_eQTL_rds,genes = c('ABCA1'),labels = labels[2])
	galactose_eQTL_plot = plot_treatment_specific_eQTL(galactose_eQTL_rds,genes = c('PRPF8'),labels = labels[3])
	shared_eQTL_plot = plot_treatment_specific_eQTL(shared_eQTL_rds,genes = c('RGR'),labels=labels[4])
	rpe_specific_eQTL = plot_rpe_specific_eQTL(rpe_specific_eQTL_rds,labels=labels[5:7],out_dir)

	topright = plot_grid(glucose_eQTL_plot,galactose_eQTL_plot,shared_eQTL_plot,nrow=3,align='v',labels = c('','',''))
	top = plot_grid(blank,topright,nrow=1,rel_widths = c(4,2.5), align = 'h',labels = c(labels[1],''))
	entire = plot_grid(top,rpe_specific_eQTL,nrow=2,align='v',labels = c('',''),rel_heights=c(4,3))
	return(entire)
}

entire = plot_entire_figure(labels=letters[1:7])
fig_fn = sprintf('%s/figure4_missingA_lc.pdf',fig_dir)
save_plot(fig_fn,entire,base_height=8.5,base_width=7)

entire = plot_entire_figure(labels=LETTERS[1:7])
fig_fn = sprintf('%s/figure4_missingA_uc.pdf',fig_dir)
save_plot(fig_fn,entire,base_height=8.5,base_width=7)

# Output underlying data:
glucose_eQTL = readRDS(glucose_eQTL_rds)
fwrite(glucose_eQTL[['ABCA1']]$data,sprintf('%s/ABCA1.txt',out_dir),sep='\t')

galactose_eQTL = readRDS(galactose_eQTL_rds)
fwrite(galactose_eQTL[['PRPF8']]$data,sprintf('%s/PRPF8.txt',out_dir),sep='\t')

shared_eQTL = readRDS(shared_eQTL_rds)
fwrite(shared_eQTL[['RGR']]$data,sprintf('%s/RGR.txt',out_dir),sep='\t')

rpe_specific_eQTL = readRDS(rpe_specific_eQTL_rds)
rpe_specific_eQTL