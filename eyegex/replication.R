library(openxlsx)
library(data.table)
library(stringr)
source('utils/genome_annotation.R')

eyegex_fn = '../data/eyegex/eyegex-supp3.xlsx'
rpe_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenesMT.txt'

mean_expression = read_mean_expression()
gene_annotation = read_gencode()
mean_expression = merge(mean_expression,gene_annotation,by=c('gene_id','gene_name'))
expressed_genes = mean_expression[mean_rpkm>0.5,]

out_dir = '../processed_data/eyegex/replication/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

eyegex = as.data.table(read.xlsx(eyegex_fn,sheet = 1, rows = 6:14865))


# Number of genes:
length(unique(eyegex$external_gene_name)) # 10463

# P-value distribution:
hist(eyegex$backward_nom_pval,breaks = 100)

# Read rpe eQTL:
rpe = fread(rpe_fn)
rpe = rpe[gene %in% expressed_genes$gene_id]

# Adding gene name:
# id2name = fread('../data/rnaseq/mean_expression.txt')[,c(1,2)]
# rpe$gene_name = id2name[match(rpe$gene,id2name$gene_id),gene_name]
# fwrite(rpe,rpe_fn,sep='\t')

rpe = rpe[glucose == 1 & galactose == 1] # selecting shared eQTL 
rpe$gene = str_split_fixed(rpe$gene,'\\.',2)[,1]
nrow(rpe) # 687

# Replicate rpe eQTL in eyegex: 
replicated = rpe$gene %in% eyegex$gene_id
sum(replicated) # 498 replicated eQTLs 
nrow(rpe[!replicated]) # 189 unreplicated eQTLs

# Write result to disk:
out_fn = paste0(out_dir, '/unreplicated_eGene.txt')
fwrite(rpe[!replicated],out_fn,sep='\t')




