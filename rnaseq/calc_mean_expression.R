library(data.table)

rpkm_fn = '../data/rnaseq/rpkm/rpe.rnaseqc_rpkm'
count_fn = '../data/rnaseq/count/merged/rpe.gene_count'
out_fn = '../data/rnaseq/mean_expression.txt'

rpkm = fread(rpkm_fn)
mean_rpkm = apply(rpkm[,3:ncol(rpkm)],1,mean)
mean_rpkm = data.table(gene_id = rpkm$Name, gene_name = rpkm$Description, mean_rpkm = mean_rpkm)


count = fread(count_fn)
mean_count = apply(count[,2:ncol(count)],1,mean)
mean_count = data.table(gene_id = count$gene_id, mean_count = mean_count)

out = merge(mean_rpkm,mean_count,by='gene_id')

fwrite(out,out_fn,sep='\t')