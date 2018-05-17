# Code to produce figure 1
# Boxiang Liu
# 2018/05/16
library(data.table)
library(manhattan)
library(cowplot)
library(ggrepel)

gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
fig_dir = '../figures/figure1/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

make_plot_data = function(out,biotype){
	plot_data = copy(out)
	plot_data[,c('chrom','pos','y'):=list(chr,start,zscore)]
	plot_data = plot_data[!chr%in%c('chrX','chrY','chrM'),]
	plot_data[,color:=ifelse(abs(zscore)>4,'red',NA)]
	plot_data[color=='red',color:=ifelse(chr%in%paste0('chr',seq(1,22,2)),'red','pink')]
	mhc_genes = plot_data[chr=='chr6'&start>28477797&stop<33448354,gene_id]
	plot_data = plot_data[!gene_id%in%mhc_genes,] # remove MHC genes
	plot_data = plot_data[type==biotype]
	plot_data[,rank := rank(-zscore)]
	plot_data[,label := ifelse(rank<=15,gene_name,NA)]
	plot_data = add_cumulative_pos(plot_data, 'hg19')
	return(plot_data)
}

plot_zscore = function(plot_data,top){
	p = manhattan(plot_data,build = 'hg19') + 
		geom_text_repel(
			data = plot_data[rank<=top],
			aes(label = gene_name),
			color='black',
			direction='x',
			segment.color = "grey",
			nudge_y=9-plot_data[rank<=top]$y,
			angle=90)+
		ylab('RPE specificity (z-score)')+
		scale_y_continuous(breaks = seq(-2,9,2),limits = c(-3,9))
	return(p)
}

read_rpe_specific_gene = function(rpe_specific_gene_fn){
	rpe_specific_gene = fread(rpe_specific_gene_fn)
	rpe_specific_gene = rpe_specific_gene[!chr%in%c('chrX','chrY','chrM'),]
	mhc_genes = rpe_specific_gene[chr=='chr6'&start>28477797&stop<33448354,gene_id]
	rpe_specific_gene = rpe_specific_gene[!gene_id%in%mhc_genes,] # remove MHC genes
	setorder(rpe_specific_gene,-zscore)
	return(rpe_specific_gene)
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


plot_gene = function(plot_data,tissue_color,tissue_abbreviation,title,top){
	median_rpkm = plot_data[,list(rpkm = median(rpkm)),by='tissue']
	setorder(median_rpkm,-rpkm)
	median_rpkm = median_rpkm[1:top,]
	plot_data = plot_data[tissue %in% median_rpkm$tissue,]
	plot_data$tissue = factor(plot_data$tissue,levels = rev(median_rpkm$tissue))
	p = ggplot(plot_data,aes(tissue,rpkm,fill = tissue)) + 
		geom_boxplot(outlier.size=-1) + 
		scale_x_discrete(breaks = names(tissue_abbreviation), labels = tissue_abbreviation) + 
		scale_fill_manual(values = tissue_color, guide = 'none') + 
		xlab('') + 
		ylab('RPKM') + 
		scale_y_log10(breaks=seq(1,10,2)) + 
		coord_flip() + 
		ggtitle(title)
	return(p)
}

# MDS plot: 
mds_fn = '../processed_data/mds/mds.GTExV7/small_mds.rds'
load(mds_fn)
blank = ggplot() + geom_blank()

# z-score Manhattan plot:
zscore_fn = '../processed_data/rpe_specific_genes.GTExV7/all_genes.txt'
zscore = fread(zscore_fn)
plot_data = make_plot_data(zscore,'protein_coding')
p = plot_zscore(plot_data,top=17)

# RPE-specific gene examples:
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)
tissue_color = read_tissue_color(gtex_tissue_color_fn)
top = 25
rpe_specific_gene_dir = '../processed_data/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
gene_name1 = 'RPE65'
fn1 = sprintf('%s/%s.txt',rpe_specific_gene_dir,gene_name1)
x1 = fread(fn1)
p1 = plot_gene(x1,tissue_color,tissue_abbreviation,title=gene_name1,top=top) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

rpe_specific_gene_dir = '../processed_data/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
gene_name2 = 'RGR'
fn2 = sprintf('%s/%s.txt',rpe_specific_gene_dir,gene_name2)
x2 = fread(fn2)
p2 = plot_gene(x2,tissue_color,tissue_abbreviation,title=gene_name2,top =top) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

left = plot_grid(blank,p,nrow=2,rel_heights = c(2,1),labels = c('A','B'))
right = plot_grid(p1,p2,nrow=2,rel_heights = c(1,1),labels = c('C','D'))
entire = plot_grid(left,right,ncol=2,rel_widths = c(3,1),labels = '')
save_plot(sprintf('%s/figure1.pdf',fig_dir),entire,base_height=9,base_width=8)

