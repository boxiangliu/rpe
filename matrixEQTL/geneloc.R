gene_annotation=read.table('../data/gtex/gencode.v19.genes.v6p.patched_contigs_genetypes.bed',header=F,stringsAsFactors=F,col.names=c('chr','start','end','strand','gene_id','type'))
subset=read.table('../processed_data/response_eQTL/residual/residual_difference.tsv',stringsAsFactors=FALSE,header=TRUE)[,1]

idx=match(subset,gene_annotation$gene_id)
out=gene_annotation[idx,c('gene_id','chr','start','end')]
out$chr=paste0('chr',out$chr)

write.table(out,'../processed_data/matrixEQTL/geneloc/geneloc.txt',sep='\t',row.names=FALSE,quote=FALSE)