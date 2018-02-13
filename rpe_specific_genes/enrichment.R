library(data.table)
fig_dir='../figures/rpe_specific_genes/enrichment/'
if (!dir.exists(fig_dir)) dir.create(fig_dir)

rpe_specific_genes_fn = '../processed_data/rpe_specific_genes/all_genes.txt'
rpe_specific_genes_qnorm_fn = '../processed_data/rpe_specific_genes/all_genes_qnorm.txt'
known_rpe_markers_fn = 'rpe_marker_gene/rpe_signature.txt'

rpe_specific_genes = fread(rpe_specific_genes_fn)[type=='protein_coding'&!is.na(zscore)]
setorder(rpe_specific_genes,-zscore)

rpe_specific_genes_qnorm = fread(rpe_specific_genes_qnorm_fn)[type=='protein_coding'&!is.na(zscore)]
setorder(rpe_specific_genes_qnorm,-zscore)

known_rpe_marker = fread(known_rpe_markers_fn)[,list(gene_name=`Gene Symbol`)]

num_genes_in_known_marker=cumsum(rpe_specific_genes$gene_name %in%known_rpe_marker$gene_name)
num_genes_in_known_marker_qnorm=cumsum(rpe_specific_genes_qnorm$gene_name %in%known_rpe_marker$gene_name)

pdf(sprintf('%s/enrichment.pdf',fig_dir))
plot(num_genes_in_known_marker,pch=16,xlab='Z-score rank',ylab='Number of genes in known marker')
points(num_genes_in_known_marker_qnorm,col='red',pch=16)
legend('topleft',legend=c('raw','quant-norm'),col=c('black','red'),pch=16)
dev.off()