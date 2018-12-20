library(data.table)
library(stringr)

glucose_eQTL_fn = '../processed_data/table1/table1/glucose_eQTL_expressed.txt'
galactose_eQTL_fn = '../processed_data/table1/table1/galactose_eQTL_expressed.txt'
glucose_eQTL = fread(glucose_eQTL_fn)
galactose_eQTL = fread(galactose_eQTL_fn)

diffexpr_fn = '/users/nsabell/rpe/diffexp/final/filt/DESeq2_diff_exp_allCov.txt'
diffexpr = read.table(diffexpr_fn)
diffexpr$gene_name = rownames(diffexpr)
diffexpr = data.table(diffexpr)

gencode = fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf', skip = 6)
gencode = gencode[V3 == 'gene']
gencode[,gene_id:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
gencode[,gene_name:=str_extract(V9,'(?<=gene_name ")(.+?)(?=";)')]



glucose_eQTL = merge(glucose_eQTL, diffexpr, by = 'gene_name')
nrow(glucose_eQTL) # 264
nrow(glucose_eQTL[padj < 0.05]) # 66
nrow(glucose_eQTL[padj < 0.05 & log2FoldChange > 0]) # 50
nrow(glucose_eQTL[padj < 0.05 & log2FoldChange < 0]) # 16


galactose_eQTL = merge(galactose_eQTL, diffexpr, by = 'gene_name')
nrow(galactose_eQTL) # 166
nrow(galactose_eQTL[padj < 0.05]) # 37
nrow(galactose_eQTL[padj < 0.05 & log2FoldChange > 0]) # 15
nrow(galactose_eQTL[padj < 0.05 & log2FoldChange < 0]) # 22