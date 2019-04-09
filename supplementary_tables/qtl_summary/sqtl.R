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
glucose_perm_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
galactose_perm_fn = '../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz'
glucose_top_intron_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/glucose/top_intron.txt'
galactose_top_intron_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/galactose/top_intron.txt'


read_fastqtl_permutation_mode = function(fn){
	fastqtl = read.table(fn,header=FALSE,stringsAsFactors=FALSE)[,c(1,6,13,16)]
	setDT(fastqtl)
	setnames(fastqtl,c('intron','snp','beta','pval'))
	return(fastqtl)
}

read_sqtl_result = function(fn){
	sqtl = fread(fn,header=TRUE)
	appendix = foreach(i = seq(nrow(sqtl)),.combine='rbind')%dopar%{
		if (str_detect(sqtl$gene_id[i],',')){
			gene_id = str_split(sqtl$gene_id[i],',')[[1]]
			data.table(gene_id = gene_id,sqtl[i,list(intron,pval,bonf,fdr)])
		} else {
			data.table()
		}
	}
	sqtl = sqtl[!str_detect(gene_id,',')]
	sqtl = rbind(sqtl,appendix)
	return(sqtl)
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

sqtl_perm_fn_list = c(glucose = glucose_perm_fn, galactose = galactose_perm_fn)
top_intron_fn_list = c(glucose = glucose_top_intron_fn,galactose = galactose_top_intron_fn)
gene_annotation = read_gene_annotation()

for (i in seq_along(sqtl_perm_fn_list)){
	sqtl_perm_fn = sqtl_perm_fn_list[i]
	top_intron_fn = top_intron_fn_list[i]
	condition = names(sqtl_perm_fn_list)[i]
	sqtl_perm = read_fastqtl_permutation_mode(sqtl_perm_fn)

	top_intron = read_sqtl_result(top_intron_fn)
	top_intron = classify_genes(top_intron, gene_annotation)
	top_intron[is.na(type),type:='unannotated']
	top_intron = select_one_annotation_per_gene(top_intron)
	
	sqtl_perm$cluster = extract_intron_cluster(sqtl_perm$intron)
	sqtl_perm$bonf = bonferroni_correction_by_intron_cluster(sqtl_perm)
	sqtl_perm = select_top_intron_per_cluster(sqtl_perm)
	sqtl_perm$fdr = p.adjust(sqtl_perm$bonf, method='fdr')

	significant_sqtl = sqtl_perm[fdr<=0.05]
	significant_sqtl$bonf = NULL
	significant_sqtl$fdr = NULL
	significant_sqtl$cluster = NULL

	# Select only protein_coding and lincRNA
	significant_sqtl = significant_sqtl[intron %in% top_intron[type %in% c('protein_coding','lincRNA'),intron]]
	out_fn = sprintf('%s/sqtl_%s_fdr0.05.txt',out_dir,condition)
	fwrite(significant_sqtl,out_fn,sep='\t')
}

