library(cowplot)
library(data.table)
library(foreach)

response_eQTL_dir = '../processed_data/QTL_landscape/response_eQTL/'
glucose_eQTL_rds = paste0(response_eQTL_dir,'glucose_specific_eQTL_plots.rds')
galactose_eQTL_rds = paste0(response_eQTL_dir,'galactose_specific_eQTL_plots.rds')
shared_eQTL_rds = paste0(response_eQTL_dir,'shared_eQTL_plots.rds')
rpe_specific_eQTL_rds = '../processed_data/specific_eQTL/specific_eGenes_v2/rpe_specific_eGenes.rds'
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
fig_dir = '../figures/figure4/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

plot_blank = function(){
	ggplot() + geom_blank()
}

plot_treatment_specific_eQTL = function(rds, genes, labels){
	eQTL = readRDS(rds)
	p1 = eQTL[[genes[1]]] + theme(strip.text = element_blank(),strip.background=element_blank(),axis.text.x=element_text(size=7)) + ylab('')
	p2 = eQTL[[genes[2]]] + theme(strip.text = element_blank(),strip.background=element_blank(),axis.text.x=element_text(size=7)) + ylab('')
	plot_grid(p1,p2,nrow = 1, align = 'h', labels = labels) 
}


read_tissue_color = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_color = gtex_tissue_color$tissue_color_hex
	names(tissue_color) = gtex_tissue_color$tissue_site_detail_id
	tissue_color = c(tissue_color,c(`RPE - glucose` = '#FF0000', `RPE - galactose` = '#00FF00'))
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
		geom_bar(stat = 'identity') + 
		geom_hline(yintercept = -log10(0.1),color='black',linetype='dashed') +
		scale_fill_manual(guide='none',values = tissue_color) + 
		scale_x_discrete(breaks = names(tissue_abbreviation), labels = tissue_abbreviation) + 
		ylab('- Log10(FDR)') + 
		xlab('') + 
		coord_flip()
}

plot_rpe_specific_eQTL = function(rpe_specific_eQTL_rds,labels = c('H','I','J')){
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
		plot_FDR(data) + ggtitle(gene_name)
	}

	col2 = plot_grid(p[[2]],p[[3]],align='v',nrow=2,labels = labels[2:3],rel_heights = c(5.5,8))
	grid_p = plot_grid(p[[1]],col2,nrow=1,labels=c(labels[1],''))
	return(grid_p)
}

blank = plot_blank()
glucose_eQTL_plot = plot_treatment_specific_eQTL(glucose_eQTL_rds,genes = c('ANO3','LAMA4'),labels = c('B','C'))
galactose_eQTL_plot = plot_treatment_specific_eQTL(galactose_eQTL_rds,genes = c('AKNA','MCM6'),labels = c('D','E'))
shared_eQTL_plot = plot_treatment_specific_eQTL(shared_eQTL_rds,genes = c('HLA-DQB1','SPATC1L'),labels=c('F','G'))
rpe_specific_eQTL = plot_rpe_specific_eQTL(rpe_specific_eQTL_rds)

topright = plot_grid(glucose_eQTL_plot,galactose_eQTL_plot,shared_eQTL_plot,nrow=3,align='v',labels = c('','',''))
top = plot_grid(blank,topright,nrow=1,align = 'h',labels = c('A',''))
entire = plot_grid(top,rpe_specific_eQTL,nrow=2,align='v',labels = c('',''),rel_heights=c(4,3))
fig_fn = sprintf('%s/fig4_missingA.pdf',fig_dir)
save_plot(fig_fn,entire,base_height=8,base_width=8)

