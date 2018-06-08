library(data.table)
source('utils/genome_annotation.R')

GEDi_fn = '../data/eye_disease/GEDi.txt'
eQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'

gene_annotation = read_gencode()
GEDi = fread(GEDi_fn)
eQTL = fread(eQTL_fn)
setnames(eQTL,'gene','gene_id')
eQTL = merge(eQTL,gene_annotation[,list(gene_id,gene_name,gene_type)],by='gene_id')
GEDi_eQTL = eQTL[gene_name %in% GEDi$gene_name]