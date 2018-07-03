library(data.table)
library(stringr)
library(leafcutter)
source('utils/genome_annotation.R')
source('utils/handling_vcf.R')

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


fastqtl_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
exon_file = '/srv/persistent/bliu2/tools/leafcutter/leafcutter/data/gencode19_exons.txt.gz'
out_dir = '../processed_data/sqtl/visualization/prepare_input/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


# Read fastqtl:
fastqtl = fread(sprintf('zcat %s',fastqtl_fn),select=c(1,6,7,13,16),col.names=c('intron','snp','dist','beta','pval'))
intron_split = fastqtl[,str_split_fixed(intron,':',4)]
fastqtl$cluster = paste(intron_split[,1],intron_split[,4],sep=':')
fastqtl[,fwer := p.adjust(pval,method='bonferroni'),by='cluster']
fastqtl[,rank := rank(fwer,ties.method='first'),by='cluster']
fastqtl[,df := .N - 1,by='cluster']
fastqtl = fastqtl[intron %in% c('chr12:56115278:56115473:clu_1202','chr12:56115278:56117670:clu_1202','chr12:56115731:56117670:clu_1202')]


# Write significance file:
significance = fastqtl[rank==1,list(cluster,intron,df,p=fwer)]
significance$gene_name = assign_gene(exon_file,significance$intron)
significance$intron = NULL
significance$status = 'Success'
significance$loglr = NA
significance[,p.adjust := p.adjust(p,method='fdr')]
setcolorder(significance,c(1,5,6,2,3,7,4))
out_fn = sprintf('%s/significance.txt',out_dir)
fwrite(significance,out_fn,sep='\t')

# Write effect size file:
effect_size = fastqtl[,list(intron,deltapsi = beta)]
effect_size$logef = NA
setcolorder(effect_size,c(1,3,2))
out_fn = sprintf('%s/effect_size.txt',out_dir)
fwrite(effect_size,out_fn,sep='\t')

# Write metadata file:
region = 'chr12:56115778'
vcf = '../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr12.all_filters.vcf.new.gz'
dosage = extract_dosage(region,vcf)
temp = dosage[,5:ncol(dosage)]
meta_data = data.table(names(temp),unlist(temp))
meta_data[,V2:=ifelse(V2==0,'CC','Ca')]
out_fn = sprintf('%s/meta_data.txt',out_dir)
fwrite(meta_data,out_fn,sep='\t',col.names=FALSE)
