library(data.table)
library(cowplot)
library(foreach)

rpe_specific_gene_fn = '../processed_data/rpe_specific_genes.GTExV7/specific_genes.RPE.txt'
rpkm_fn = '../processed_data/mds/preprocess.GTExV7/combined_logxp2.rpkm'
col_data_fn='../processed_data/mds/preprocess.GTExV7/combined.col'
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
fig_dir = '../figures/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
out_dir = '../processed_data/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_rpe_specific_gene = function(rpe_specific_gene_fn){
	rpe_specific_gene = fread(rpe_specific_gene_fn)
	rpe_specific_gene = rpe_specific_gene[!chr%in%c('chrX','chrY','chrM'),]
	mhc_genes = rpe_specific_gene[chr=='chr6'&start>28477797&stop<33448354,gene_id]
	rpe_specific_gene = rpe_specific_gene[!gene_id%in%mhc_genes,] # remove MHC genes
	setorder(rpe_specific_gene,-zscore)
	return(rpe_specific_gene)
}

read_rpkm = function(rpkm_fn){
	rpkm = fread(rpkm_fn)
	return(rpkm)
}

read_col_data = function(col_data_fn){
	col_data = fread(col_data_fn)
	return(col_data)
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

get_tissue_to_tissue_id = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	tissue_to_tissue_id = gtex_tissue_color$tissue_site_detail_id
	names(tissue_to_tissue_id) = gtex_tissue_color$tissue_site_detail
	tissue_to_tissue_id = c(tissue_to_tissue_id,`RPE (glu)`="RPE - glucose",`RPE (gal)`="RPE - galactose")
	return(tissue_to_tissue_id)
}

make_plot_data = function(rpkm,col_data,gene_id){
	gene_rpkm = data.frame(rpkm = t(rpkm[Name == gene_id,2:ncol(rpkm)]))
	gene_rpkm$sample = rownames(gene_rpkm)
	setDT(gene_rpkm)
	plot_data = merge(gene_rpkm,col_data,by='sample')
	median_rpkm = plot_data[,list(median = median(rpkm)),by='tissue']
	setorder(median_rpkm,'median')
	plot_data$tissue = factor(plot_data$tissue,levels = median_rpkm$tissue)
	return(plot_data)
}

plot_gene = function(plot_data,tissue_color,tissue_abbreviation,title){
	p = ggplot(plot_data,aes(tissue,rpkm,fill = tissue)) + 
		geom_boxplot(outlier.size=-1) + 
		scale_x_discrete(breaks = names(tissue_abbreviation), labels = tissue_abbreviation) + 
		scale_fill_manual(values = tissue_color, guide = 'none') + 
		xlab('') + 
		ylab(expression(paste(log[10],'(RPKM+2)'))) + 
		scale_y_log10() + 
		coord_flip() + 
		ggtitle(title)
	return(p)
}

make_plot_data_2 = function(rpe_specific_proteinCodingGene){
	plot_data_2 = foreach (i = 1:10,.combine='rbind')%do%{
		gene_id = rpe_specific_proteinCodingGene$gene_id[i]
		gene_name = rpe_specific_proteinCodingGene$gene_name[i]
		plot_data = make_plot_data(rpkm,col_data,gene_id = gene_id)

		x = plot_data[, list(rpkm = mean(rpkm)), by = 'tissue']
		x[,rank:=rank(rpkm)]
		x$gene = gene_name
		return(x)
	}
	return(plot_data_2)
}

plot_top_rpe_specific_genes = function(plot_data_2){
	p2 = ggplot(plot_data_2,aes(rank,rpkm,color=tissue))+
		geom_point()+
		facet_wrap(~gene,nrow=2,scales='free_y')+
		scale_color_manual(values=tissue_color,labels=tissue_abbreviation,name=NULL)+
		scale_y_log10(breaks=seq(1,10,2))+
		theme_bw()+
		theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),panel.grid=element_blank(),strip.background = element_rect(color='black',fill='white'),legend.position='bottom') + 
		xlab('Tissues') + 
		ylab('Mean RPKM') 
	return(p2)
}

# Read data:
rpe_specific_gene = read_rpe_specific_gene(rpe_specific_gene_fn)
rpkm = read_rpkm(rpkm_fn)
col_data = read_col_data(col_data_fn)
tissue_color = read_tissue_color(gtex_tissue_color_fn)
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)
tissue_to_tissue_id = get_tissue_to_tissue_id(gtex_tissue_color_fn)
col_data$tissue = tissue_to_tissue_id[col_data$tissue]
rpe_specific_proteinCodingGene = rpe_specific_gene[type=='protein_coding']

# Make main figure:
for (i in 1:20){
	gene_id = rpe_specific_proteinCodingGene$gene_id[i]
	gene_name = rpe_specific_proteinCodingGene$gene_name[i]

	plot_data = make_plot_data(rpkm,col_data,gene_id = gene_id)
	p = plot_gene(plot_data,tissue_color,tissue_abbreviation,title=gene_name)

	fwrite(plot_data,sprintf('%s/%s.txt',out_dir,gene_name),sep='\t')
	save_plot(sprintf('%s/%s.pdf',fig_dir,gene_name),p,base_height=8,base_width=4)
}


# Make supplementary figure:
plot_data_2 = make_plot_data_2(rpe_specific_proteinCodingGene)
p2 = plot_top_rpe_specific_genes(plot_data_2)
save_plot(sprintf('%s/top_10_rpe_specific_genes.pdf',fig_dir),p2,base_width=8.5,base_height=7)