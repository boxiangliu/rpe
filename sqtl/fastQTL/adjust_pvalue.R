library(data.table)
library(stringr)
library(leafcutter)

glucose_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
galactose_fn = '../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz'
gtf_fn = '../data/reference/gencode.v19.annotation.gtf'
out_dir = '../processed_data/sqtl/fastQTL/adjust_pvalue/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
exon_file = '../processed_data/sqtl/fastQTL/adjust_pvalue/gencode19_exons.gene_id.txt.gz'

gtf_to_exon = function(gft_fn){
	gtf = read.table(gtf_fn,stringsAsFactors=FALSE,header=FALSE,sep='\t',quote="")[,c(1,3,4,5,7,9)]
	setDT(gtf)
	setnames(gtf,c('chr','feature','start','end','strand','annotation'))
	exon = gtf[feature == 'exon']
	if (!str_detect(exon$chr[1],'chr')) {exon[,chr := paste0('chr',chr)]}
	exon[,gene_name := str_extract(annotation,'(?<=gene_id\\s\\")(ENSG.+?)(?=\\";)')]
	exon$feature = NULL
	exon$annotation = NULL
	return(exon)
}

read_fastqtl_permutation_mode = function(fn){
	fastqtl = read.table(sqtl_fn,header=FALSE,stringsAsFactors=FALSE)[,c(1,16)]
	setDT(fastqtl)
	setnames(fastqtl,c('intron','pval'))
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

save_top_intron = function(top_intron,fn){
	top_intron$cluster = NULL
	fwrite(top_intron,fn,sep='\t')
}

exon = gtf_to_exon(gtf_fn)
write.table(exon,gzfile(exon_file),col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')

temp = c(glucose = glucose_fn, galactose = galactose_fn)

for (i in seq_along(temp)){
	sqtl_fn = temp[i]
	sqtl = read_fastqtl_permutation_mode(sqtl_fn)
	sqtl$gene_id=assign_gene(exon_file,sqtl$intron)

	sqtl$cluster = extract_intron_cluster(sqtl$intron)
	sqtl$bonf = bonferroni_correction_by_intron_cluster(sqtl[,list(cluster,pval)])

	top_intron = select_top_intron_per_cluster(sqtl)
	top_intron$fdr = p.adjust(top_intron$bonf, method='fdr')

	out_dir_2 = sprintf('%s/%s/',out_dir,names(temp)[i])
	if (!dir.exists(out_dir_2)) {dir.create(out_dir_2)}
	out_fn = paste0(out_dir_2,'/top_intron.txt')
	save_top_intron(top_intron,out_fn)
}

