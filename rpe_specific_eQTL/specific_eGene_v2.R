library(data.table)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(10)
library(cowplot)

#-----------#
# Variables #
#-----------#
treeQTL_MT_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenesMT.txt'
GTEx_v7_dir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/'
glucose_eGenes_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenes_glucose.txt'
galactose_eGenes_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenes_galactose.txt'
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
fig_dir = '../figures/rpe_specific_eQTL/specific_eGenes_v2/'
out_dir = '../processed_data/rpe_specific_eQTL/specific_eGenes_v2/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

#-----------#
# Functions #
#-----------#
read_TreeQTL_MT = function(fn){
	treeQTL_MT = fread(fn)
	return(treeQTL_MT)
}

read_TreeQTL = function(fn){
	treeQTL = fread(fn)
	return(treeQTL)
}

read_GTEx = function(fn){
	tissue = basename(fn)
	tissue = str_replace(tissue,'.v7.egenes.txt.gz','')
	command = sprintf('gunzip -c %s',fn)
	GTEx = fread(cmd = command)
	GTEx = GTEx[,list(gene_id,gene_name,variant_id,rs_id_dbSNP147_GRCh37p13,pval_beta,slope,slope_se,qval)]
	GTEx$tissue = tissue
	return(GTEx)
}

screen_GTEx = function(gene,GTEx,QVAL){
	qval = GTEx[gene_id == gene, list(qval,tissue,slope,slope_se,gene_name)]
	if (nrow(qval)==0){
		return(NULL)
	}
	gene_name = qval$gene_name
	not_eQTL = all(qval$qval > QVAL)
	min_qval = min(qval$qval)
	tissue = qval[which.min(qval),tissue]
	slope = qval[which.min(qval),slope]
	slope_se = qval[which.min(qval),slope_se]
	GTEx_eQTL = data.table(gene_id = gene, gene_name, eQTL = !not_eQTL,slope,slope_se,min_qval,tissue)
	return(unique(GTEx_eQTL))
}

read_tissue_color = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_color = gtex_tissue_color$tissue_color_hex
	names(tissue_color) = gtex_tissue_color$tissue_site_detail_id
	tissue_color = c(tissue_color,c(`RPE - glucose` = '#F87660', `RPE - galactose` = '#619CFF'))
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
#------#
# Main #
#------#
# Read data: 
treeQTL_MT = read_TreeQTL_MT(treeQTL_MT_fn)
shared_eQTL = treeQTL_MT[glucose == 1 & galactose == 1]

fn_list = list.files(GTEx_v7_dir,pattern='egenes',full.names=TRUE)

GTEx = foreach(fn = fn_list,.combine = rbind)%dopar%{
	read_GTEx(fn)
}

# Find RPE-specific eQTL (GTEx FDR > 0.1):
QVAL = 0.1
GTEx_eQTL = foreach(gene = shared_eQTL[,gene],.combine = rbind)%dopar%{
	screen_GTEx(gene,GTEx,QVAL)
}


# Make plot data:
glucose_eGenes = read_TreeQTL(glucose_eGenes_fn)
galactose_eGenes = read_TreeQTL(galactose_eGenes_fn)

tissue_color = read_tissue_color(gtex_tissue_color_fn)
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)

# Make plots:
gene_list = unique(GTEx_eQTL[eQTL == FALSE,gene_id])
p = foreach(gene = gene_list)%do%{
	data = make_plot_data(gene, GTEx, glucose_eGenes, galactose_eGenes)
	gene_name = GTEx_eQTL[gene_id==gene,unique(gene_name)]
	plot_FDR(data) + ggtitle(gene_name)
}
names(p) = gene_list
out_fn = sprintf('%s/rpe_specific_eGenes.rds',out_dir)
saveRDS(list(GTEx_eQTL=GTEx_eQTL,GTEx=GTEx,glucose_eGenes=glucose_eGenes,galactose_eGenes=galactose_eGenes,
	tissue_color=tissue_color,tissue_abbreviation=tissue_abbreviation),out_fn)

col2 = plot_grid(p[[2]],p[[3]],align='v',nrow=2,labels = c('B','C'),rel_heights = c(5.5,8))
grid_p = plot_grid(p[[1]],col2,nrow=1,labels=c('A',''))
fig_fn = sprintf('%s/rpe_specific_eGenes.pdf',fig_dir)
save_plot(fig_fn,grid_p,base_height=6,base_width=6)

# Save RPE-specific eQTL:
out_fn = sprintf('%s/rpe_specific_eGenes.txt',out_dir)
fwrite(GTEx_eQTL[eQTL == FALSE],out_fn,sep='\t')