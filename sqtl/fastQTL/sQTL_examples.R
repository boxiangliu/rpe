library(data.table)
source('utils/genome_annotation.R')

glucose_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/glucose/top_intron.txt'
galactose_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/galactose/top_intron.txt'

gene_annotation = read_gencode()
glucose = fread(glucose_fn)
glucose = merge(glucose,gene_annotation[,list(gene_id,gene_name)],by='gene_id')
glucose = glucose[!is.na(glucose$fdr)]
setorder(glucose,fdr)

galactose = fread(galactose_fn)
galactose = merge(galactose,gene_annotation[,list(gene_id,gene_name)],by='gene_id')
galactose = galactose[!is.na(galactose$fdr)]
setorder(galactose,fdr)
galactose