library(data.table)
library(manhattan)
library(cowplot)
library(ggrepel)
devtools::install_github('boxiangliu/manhattan')
source('utils/genome_annotation.R')


deseq_fn = '/users/nsabell/rpe/diffexp/final/filt/DESeq2_diff_exp_allCov.txt'
gene_annotation = read_gencode(gencode_fn)
mean_expression = read_mean_expression()
expressed_gene = mean_expression[mean_rpkm > 0.5, gene_id]
fig_dir = '../figures/diff_expression/plot_manhattan/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_deseq = function(deseq_fn){
	results = read.table(deseq_fn, header = T, na.strings = "NA")
	results = results[!is.na(results$padj),]
	results$gene_name = rownames(results)
	setDT(results)
	return(results)
}

annotate_deseq = function(results,expressed_gene){
	setDT(results)
	results_annotated = merge(results, gene_annotation, by = "gene_name")
	results_annotated = results_annotated[gene_id %in% expressed_gene]
	return(results_annotated)
}


make_plot_data = function(results_annotated,cutoff=1e-5,limit=1e-16){
	manhattan_data = results_annotated[,list(chrom=chr,pos=start,padj,gene_id,gene_name,log2FoldChange)]
	manhattan_data = manhattan_data[chrom %in% paste0('chr',1:22),]
	manhattan_data[,rank:=rank(padj,ties.method='first')]
	manhattan_data[,shape:=16]
	manhattan_data[,color:=ifelse(padj<cutoff&log2FoldChange>0&chrom%in%paste0('chr',seq(1,21,2)),'red',NA)]
	manhattan_data[,color:=ifelse(padj<cutoff&log2FoldChange>0&chrom%in%paste0('chr',seq(2,22,2)),'pink',color)]
	manhattan_data[,color:=ifelse(padj<cutoff&log2FoldChange<0&chrom%in%paste0('chr',seq(1,21,2)),'green',color)]
	manhattan_data[,color:=ifelse(padj<cutoff&log2FoldChange<0&chrom%in%paste0('chr',seq(2,22,2)),'greenyellow',color)]
	manhattan_data[,fill:=color]
	manhattan_data[,y := -log10(padj) * sign(log2FoldChange)]
	setorder(manhattan_data,chrom,pos)
	manhattan_data = add_cumulative_pos(manhattan_data, build = "hg19")
	return(manhattan_data)
}

plot_manhattan = function(plot_data,top=20){
	color_map = c('black','grey','red','pink','green','greenyellow')
	names(color_map) = color_map
	text_data = plot_data[which(plot_data$rank<=top),]
	set.seed(142)

	p = manhattan(plot_data) + 
		scale_color_manual(name='',values = color_map,breaks = c('red','green'),labels=c('Glucose upregulation','Galactose upregulation')) + 
		theme(legend.position = c(0.99,1.05),legend.justification=c('right','top'),legend.text=element_text(size=10)) +
		geom_text_repel(
			data = text_data,
			aes(label = gene_name),
			color='black',
			segment.color = "grey",
			box.padding = 0.2,
			angle=0,
			size=3) + 
		scale_y_continuous(name=expression(-log[10]*"(FDR) x sign("*beta*')'))
	return(p)
}

results = read_deseq(deseq_fn)
results_annotated = annotate_deseq(results,expressed_gene)
nrow(results_annotated[padj<1e-3,]) # 837
results_annotated_pc = results_annotated[gene_type=='protein_coding']
plot_data = make_plot_data(results_annotated_pc,cutoff = 1e-3)
p = plot_manhattan(plot_data,top=20)
fig_fn = sprintf('%s/deseq_manhattan.pdf',fig_dir)
ggsave(fig_fn,p,height=3,width=8)


