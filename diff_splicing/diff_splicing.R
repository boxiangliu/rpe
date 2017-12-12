# Perform differential splicing analysis with leafcutter
# Boxiang Liu
# 2017-12-11

library(leafcutter)
library(foreach)
library(doMC)
registerDoMC(10)
library(gap)
library(stringr)

count_fn='../data/rnaseq/leafcutter/both/cluster/diff_splice_perind_numers.counts.gz'
exon_file='/srv/persistent/bliu2/tools/leafcutter/leafcutter/data/gencode19_exons.txt.gz'
sample_table_fn='../processed_data/diff_splicing/sample_file/sample_table.txt'
out_prefix='../processed_data/diff_splicing/diff_splicing/gal_vs_glu'
fig_dir='../figures/diff_splicing/diff_splicing/'

if (!dir.exists(dirname(out_prefix))){dir.create(dirname(out_prefix),recursive=TRUE)}
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

read_x=function(sample_table_fn){
	predictor=read.table(sample_table_fn)
	x=as.integer(predictor$V2=='galactose')
}

read_confounders=function(sample_table_fn){
	predictor=read.table(sample_table_fn)
	confounders=predictor[,c('V3','V4','V5')]
	confounders$V5=(confounders$V5==F)
	confounders=as.matrix(confounders)
	return(confounders)
}

assign_gene=function(cluster_table,exon_file,counts){
	exons_table     = read.table(exon_file, header=T, stringsAsFactors = F)
	intron_meta     = get_intron_meta(rownames(counts))
	exons_table$chr = add_chr(exons_table$chr)
	intron_meta$chr = add_chr(intron_meta$chr)
	clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
	cluster_table   = merge(cluster_table, clu_gene_map, by.x="cluster", by.y="clu", all.x=TRUE)
	return(cluster_table)
}

make_QC_plots=function(cluster_table,effect_size_table,fig_dir){
	pdf(sprintf('%s/QC.pdf',fig_dir))
	hist(cluster_table$p,breaks=100)
	qqunif(cluster_table$p)
	hist(effect_size_table$logef)
	dev.off()
}


make_sashimi_plot=function(cluster_table,counts,x,num,fig_dir){
	cluster_table_ordered=cluster_table[order(cluster_table$p),]
	clusters=str_replace(rownames(counts),':[0-9]+:[0-9]+:',':')
	
	pdf(sprintf('%s/sashimi_plot.pdf',fig_dir))
	for (i in 1:num){
		y=t(counts[clusters%in%cluster_table_ordered$cluster[i],])
		gene_name=cluster_table_ordered$genes[i]
		
		make_differential_splicing_plot(y,x,
		  exons_table = exons_table, len = 500, 
		  length_transform = function(g) log(g+1), 
		  main_title = gene_name, summary_func = colMeans,
		  legend_title = "Mean count")
	}
	dev.off()
}

counts=read.table(count_fn,check.names=FALSE)
x=read_x(sample_table_fn)
# confounders=read_confounders(sample_table_fn)
exons_table = read.table(exon_file, header=T, stringsAsFactors = F)

results=differential_splicing(counts,x)

cluster_table = cluster_results_table(results)
cluster_table$cluster = add_chr(cluster_table$cluster)
cluster_table=assign_gene(cluster_table,exon_file,counts)
write.table(cluster_table, paste0(out_prefix,"_cluster_significance.txt"), quote=F, sep="\t", row.names = F)

effect_size_table = leaf_cutter_effect_sizes(results)
setnames(effect_size_table,c(''),)
colnames(effect_size_table)[3:4]=c('glu','gal')
write.table(effect_size_table,paste0(out_prefix,"_effect_sizes.txt"), quote=F, col.names = T, row.names = F, sep="\t")

# sum(cluster_table$p.adjust<0.05,na.rm=TRUE) # 21

make_QC_plots(cluster_table,effect_size_table,fig_dir)
make_sashimi_plot(cluster_table,counts,x,10,fig_dir)