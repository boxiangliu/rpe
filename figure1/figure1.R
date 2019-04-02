# Code to produce figure 1
# Boxiang Liu
# 2018/05/16
library(data.table)
library(manhattan)
library(cowplot)
library(ggrepel)
source('utils/gtex_tissue_info.R')
library(openxlsx)

gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
highlight_gene_fn = 'figure1/Fig 1b genes.xlsx'
fig_dir = '../figures/figure1/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
out_dir = '../processed_data/figure1/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

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

plot_zscore = function(plot_data,highlight=NULL,top=16){
	plot_data = add_shape(plot_data)
	plot_data = add_fill(plot_data)
	if (!is.null(highlight)){
		p = manhattan(plot_data,build = 'hg19') + 
			geom_text_repel(
				data = plot_data[gene_name %in% highlight],
				aes(label = gene_name),
				color='black',
				direction='x',
				segment.color = "grey",
				nudge_y=10-plot_data[gene_name %in% highlight]$y,
				angle=90)+
			ylab('RPE specificity (z-score)')+
			scale_y_continuous(breaks = seq(-2,10,2),limits = c(-3,10))
	} else {
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
	}
	return(p)
}

read_highlight_gene = function(fn){
	x = read.xlsx(fn,sheet = 1)
	setDT(x)
	chr1 = c('RPE65','PLD5','CA14','ABCA4')
	chr2 = c('TTLL4')
	chr3 = c('SLC6A20')
	chr4 = c('LRAT')
	chr5 = c('SLC45A2')
	chr6 = c('SLC2A12')
	chr7 = c('SLC4A2','STRIP2')
	chr8 = c()
	chr9 = c('TYRP1','TRPM3')
	chr10 = c('RGR')
	chr11 = c('BEST1','TYR','MFRP')
	chr12 = c('RDH5')
	chr13 = c('DCT')
	chr14 = c('OTX2')
	chr15 = c('RLBP1')
	chr16 = c('BCMO1')
	chr17 = c()
	chr18 = c('RAX')
	chr19 = c('CRX')
	chr20 = c()
	chr21 = c()
	chr22 = c('SLC16A8')
	select = c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
	x = x[gene_name %in% select]
	return(x)
}

read_rpe_specific_gene = function(rpe_specific_gene_fn){
	rpe_specific_gene = fread(rpe_specific_gene_fn)
	rpe_specific_gene = rpe_specific_gene[!chr%in%c('chrX','chrY','chrM'),]
	mhc_genes = rpe_specific_gene[chr=='chr6'&start>28477797&stop<33448354,gene_id]
	rpe_specific_gene = rpe_specific_gene[!gene_id%in%mhc_genes,] # remove MHC genes
	setorder(rpe_specific_gene,-zscore)
	return(rpe_specific_gene)
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
		ylab(expression(paste(log[2],'(RPKM+2)'))) + 
		scale_y_log10(breaks=seq(1,10,2)) + 
		coord_flip() + 
		ggtitle(title)
	return(p)
}


# Read tissue color and abbreviation:
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)
tissue_color = read_tissue_color(gtex_tissue_color_fn)

# MDS plot: 
mds_fn = '../processed_data/mds/mds.GTExV7/mds_plot_data.rda'
result = readRDS(mds_fn)
mds = result[[1]]
fwrite(mds, '../processed_data/figure1/mds.txt', sep = '\t')
mds_centroid = result[[2]]
p1=ggplot(mds,aes(x=x,y=y,color=tissue))+
	geom_point()+scale_color_manual(values=tissue_color,guide='none')+
	geom_text(data=mds_centroid,aes(x=x_centroid,
		y=y_centroid,label=label),color='black')+
	theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=12))


# z-score Manhattan plot:
zscore_fn = '../processed_data/rpe_specific_genes/rpe_specific_genes.GTExV7/expressed_genes.RPE.txt'
zscore = fread(zscore_fn)
plot_data = make_plot_data(zscore,'protein_coding')
highlight_gene = read_highlight_gene(highlight_gene_fn)
set.seed(102)
p2 = plot_zscore(plot_data,highlight=highlight_gene$gene_name) + theme(axis.title=element_text(size=12))


# RPE-specific gene examples:
top = 25
rpe_specific_gene_dir = '../processed_data/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
gene_name1 = 'RPE65'
fn1 = sprintf('%s/%s.txt',rpe_specific_gene_dir,gene_name1)
x1 = fread(fn1)
p3 = plot_gene(x1,tissue_color,tissue_abbreviation,title=gene_name1,top=top) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_text(size=12))


rpe_specific_gene_dir = '../processed_data/rpe_specific_genes/plot_rpe_specific_gene.GTExV7/'
gene_name2 = 'RGR'
fn2 = sprintf('%s/%s.txt',rpe_specific_gene_dir,gene_name2)
x2 = fread(fn2)
p4 = plot_gene(x2,tissue_color,tissue_abbreviation,title=gene_name2,top =top) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title=element_text(size=12))

left = plot_grid(p1,p2,nrow=2,rel_heights = c(2,1),labels = c('A','B'))
right = plot_grid(p3,p4,nrow=2,rel_heights = c(1,1),labels = c('C','D'))
entire = plot_grid(left,right,ncol=2,rel_widths = c(3,1),labels = '')
save_plot(sprintf('%s/figure1_uc.pdf',fig_dir),entire,base_height=8.5,base_width=8)


left = plot_grid(p1,p2,nrow=2,rel_heights = c(2,1),labels = c('a','b'))
right = plot_grid(p3,p4,nrow=2,rel_heights = c(1,1),labels = c('c','d'))
entire = plot_grid(left,right,ncol=2,rel_widths = c(3,1),labels = '')
save_plot(sprintf('%s/figure1_lc.pdf',fig_dir),entire,base_height=8.5,base_width=8)

