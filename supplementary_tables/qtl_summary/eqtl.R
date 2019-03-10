library(data.table)
library(stringr)
library(foreach)
library(doMC)
registerDoMC(15)
source('utils/genome_annotation.R')

treeQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenesMT.txt'
glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
out_dir = '../processed_data/supplementary_tables/qtl_summary/eqtl/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

treeQTL = fread(treeQTL_fn)
glucose = treeQTL[glucose==1,gene]
galactose = treeQTL[galactose==1,gene]

select_top_eQTL = function(fn){
	eqtl = fread(fn,select = c(1,2,11,12),col.names=c('fid','snp','chisq','pi'))
	eqtl[,rank := rank(-chisq,ties.method='first')]
	eqtl = eqtl[rank==1]
	eqtl$rank = NULL
	split_fid = str_split_fixed(eqtl$fid,'_',2)
	eqtl$gene_id = split_fid[,1]
	eqtl$gene_name = split_fid[,2]
	eqtl[,pval:=pchisq(chisq,df=1,lower.tail=FALSE)]
	eqtl = eqtl[,list(gene_id,gene_name,snp,pi,pval)]
	return(eqtl)
}


mean_expression = read_mean_expression()
gene_annotation = read_gencode()
mean_expression = merge(mean_expression,gene_annotation,by=c('gene_id','gene_name'))
expressed_genes = mean_expression[mean_rpkm>0.5,]

glucose_eqtl = foreach(i = seq_along(glucose),.combine='rbind')%dopar%{
	print(i)
	fn = list.files(glucose_dir,glucose[i],recursive=TRUE,full.names=TRUE)
	eqtl = select_top_eQTL(fn)
	return(eqtl)
}

glucose_eqtl = glucose_eqtl[gene_id %in% expressed_genes$gene_id]
out_fn = sprintf('%s/glucose_eqtl_fdr0.05.txt',out_dir)
fwrite(glucose_eqtl,out_fn,sep='\t')

galactose_eqtl = foreach(i = seq_along(galactose),.combine='rbind')%dopar%{
	print(i)
	fn = list.files(galactose_dir,galactose[i],recursive=TRUE,full.names=TRUE)
	eqtl = select_top_eQTL(fn)
	return(eqtl)
}

galactose_eqtl = galactose_eqtl[gene_id %in% expressed_genes$gene_id]
out_fn = sprintf('%s/galactose_eqtl_fdr0.05.txt',out_dir)
fwrite(galactose_eqtl,out_fn,sep='\t')