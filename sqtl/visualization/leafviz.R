library(dplyr)
library(ggplot2)
library(leafcutter)
library(reshape2)
library(gridExtra)
library(intervals) # needed for pretty strand arrow placement
library(foreach)
library(grid)
library(gtable)
library(ggrepel)
library(cowplot)
source('sqtl/visualization/make_cluster_plot.R')

fig_dir = '../figures/sqtl/visualization/leafviz/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
out_dir = '../processed_data/sqtl/visualization/leafviz/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
sqtl_fn = '/srv/persistent/bliu2/rpe/processed_data/sqtl/visualization/prepare_results/sqtl.RData'
load(sqtl_fn)

filter_intron_table = function(introns, clu, toSave=FALSE){
	d = dplyr::filter(introns, clusterID == clu) %>% 
		dplyr::select( -clusterID, -gene, -ensemblID, -transcripts) %>% 
		arrange( desc(abs(deltapsi)))
	if( !toSave ){
		d = rename(d, "ΔPSI" = deltapsi )
	}else{
		d = rename(d, "dPSI" = deltapsi ) # fudge as grid arrange doesn't like greek letters
	}
	row.names(d) = letters[1:nrow(d)] # letters is just a:z
	return(d)
}

getGeneLength = function(gene_name){
	# gets length of gene in nucleotides and decides on a pixel length for the gene plot
	# RBFOX1 is 1.7Mbp - scale to 5000px
	# most genes are < 100kb
	exons = exons_table[ exons_table$gene_name == gene_name, ]
	geneStart = min(exons$start)
	geneEnd = max(exons$end)
	geneLength = geneEnd - geneStart
	#print(geneLength)
	if( geneLength >1E6){
		pixels = 5000 # scales RBFOX1 to 5000px
	}
	if( geneLength > 5e5 & geneLength < 1e6){
		pixels = 3000
	}
	if( geneLength > 1.5e5 & geneLength <= 5e5){
		pixels = 2000
	}
	if( geneLength <= 1.5e5){
		pixels = "auto"
	}
	#print(pixels)
	return(pixels)
}

get_data = function(clusters,clusterID){
	sel = which(clusters$clusterID == clusterID)
	gene = clusters[ sel, ]$gene
	gene = gsub("<.*?>", "", gene) # strip out html italic tags
	width = getGeneLength(gene)
	clusterID = clusters[ sel, ]$clusterID
	coord = clusters[ sel, ]$coord
	return(list(gene = gene, width = width, cluster = clusterID, coord = coord))
}

mydata = get_data(clusters,'clu_1202')

plotTitle = c(mydata$gene, as.character(mydata$cluster))
introns$verdict[introns$verdict == 'novel annotated pair'] = 'annotated'
introns$verdict[!(introns$start == 56115278 & introns$end == 56117670)] = 'other' 

plots = make_cluster_plot(mydata$cluster,
	main_title = NULL,
	meta = meta,
	cluster_ids = cluster_ids,
	exons_table = exons_table,
	counts = counts,
	introns = introns)
plots[[1]] = plots[[1]] 
plots[[2]] = plots[[2]] + guides(color=FALSE) 
p = plot_grid(plots[[1]],plots[[2]],nrow=2,labels='')
fig_fn = sprintf('%s/rdh5.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)
out_fn = sprintf('%s/rdh5.rda',out_dir)
saveRDS(plots,out_fn)
