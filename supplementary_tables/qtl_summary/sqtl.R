library(data.table)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(15)
library(leafcutter)
source('utils/summary_table.R')
source('utils/genome_annotation.R')

exon_file = '../processed_data/sqtl/fastQTL/adjust_pvalue/gencode19_exons.gene_id.txt.gz'
out_dir = '../processed_data/supplementary_tables/qtl_summary/sqtl/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
glucose_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
galactose_fn = '../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz'


read_fastqtl_permutation_mode = function(fn){
	fastqtl = read.table(fn,header=FALSE,stringsAsFactors=FALSE)[,c(1,6,13,16)]
	setDT(fastqtl)
	setnames(fastqtl,c('intron','snp','beta','pval'))
	return(fastqtl)
}

extract_intron_cluster = function(intron){
	intron_cluster = str_split_fixed(intron,':',4)[,4]
	return(intron_cluster)
}

bonferroni_correction_by_intron_cluster = function(x){
	x = copy(x)
	x[,n := .N, by = 'cluster']
	x[,bonf := pval*n]
	x[,bonf := ifelse(bonf > 1, 1, bonf)]
	return(x$bonf)
}

select_top_intron_per_cluster = function(sqtl){
	sqtl[, bonf_rank := rank(bonf,ties.method = 'first'), by = 'cluster']
	top_intron = sqtl[bonf_rank == 1]
	top_intron$bonf_rank = NULL
	return(top_intron)
}

assign_gene=function(exon_file,clusters){
	exons_table     = read.table(exon_file, header=T, stringsAsFactors = F)
	intron_meta     = get_intron_meta(clusters)
	exons_table$chr = add_chr(exons_table$chr)
	intron_meta$chr = add_chr(intron_meta$chr)
	clu_gene_map    = map_clusters_to_genes(intron_meta, exons_table)
	clu_map = data.frame(clusters,clu=gsub(':[0-9]+:[0-9]+','',clusters))
	clu_gene_map = merge(clu_gene_map,clu_map,by='clu')
	map = clu_gene_map$genes
	names(map) = clu_gene_map$clusters
	return(map[clusters])
}

read_gene_annotation = function(){
	gene_annotation = read_gencode(gencode_fn)
	gene_annotation[,type:=NULL]
	setnames(gene_annotation,'gene_type','type')
	return(gene_annotation)
}

select_one_annotation_per_gene = function(sqtl){
	sqtl[,type:=factor(type,levels=c('protein_coding','lincRNA','pseudogene','polymorphic_pseudogene','processed_transcript','sense_intronic','sense_overlapping','antisense','3prime_overlapping_ncrna','unannotated'))]
	sqtl[,rank:=rank(type,ties.method='first'),by='intron']
	sqtl = sqtl[rank==1]
	return(sqtl)
}

fn_list = c(glucose = glucose_fn, galactose = galactose_fn)

for (i in seq_along(fn_list)){
	fn = fn_list[i]
	condition = names(fn_list)[i]
	sqtl = read_fastqtl_permutation_mode(fn)
	sqtl$cluster = extract_intron_cluster(sqtl$intron)
	sqtl$bonf = bonferroni_correction_by_intron_cluster(sqtl)
	sqtl = select_top_intron_per_cluster(sqtl)
	sqtl$fdr = p.adjust(sqtl$bonf, method='fdr')
	significant_sqtl = sqtl[fdr<=0.05]
	significant_sqtl$bonf = NULL
	significant_sqtl$fdr = NULL
	significant_sqtl$cluster = NULL
	out_fn = sprintf('%s/sqtl_%s_fdr0.05.txt',out_dir,condition)
	fwrite(significant_sqtl,out_fn,sep='\t')
}

