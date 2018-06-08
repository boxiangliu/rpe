library(data.table)
source('utils/genome_annotation.R')

treeQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
treeQTL = fread(treeQTL_fn)
out_dir = '../processed_data/table1/table1/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

count_features = function(expressed_genes){
	tested = table(expressed_genes[,gene_type])
	protein_coding = tested['protein_coding']
	lincRNA = tested['lincRNA']
	total = protein_coding + lincRNA
	x = data.table(
		type = c('protein_coding','lincRNA','total'),
		count = c(protein_coding,lincRNA,total))
	return(x)
}

mean_expression = read_mean_expression()
gene_annotation = read_gencode()
mean_expression = merge(mean_expression,gene_annotation,by=c('gene_id','gene_name'))
expressed_genes = mean_expression[mean_rpkm>0.5,]

tested_count = count_features(expressed_genes)
setnames(tested_count,'count','tested')

glucose_eQTL = treeQTL[glucose==1&galactose==0,gene]
glucose_count = count_features(expressed_genes[gene_id%in%glucose_eQTL])
setnames(glucose_count,'count','glucose')

shared_eQTL = treeQTL[glucose==1&galactose==1,gene]
shared_count = count_features(expressed_genes[gene_id%in%shared_eQTL])
setnames(shared_count,'count','shared')

galactose_eQTL = treeQTL[glucose==0&galactose==1,gene]
galactose_count = count_features(expressed_genes[gene_id%in%galactose_eQTL])
setnames(galactose_count,'count','galactose')

merged = merge(tested_count,glucose_count,by='type')
merged = merge(merged,galactose_count,by='type')
merged = merge(merged,shared_count,by='type')

merged[,glucose_pct := sprintf('%0.2f%%',glucose/tested*100)]
merged[,galactose_pct := sprintf('%0.2f%%',galactose/tested*100)]
merged[,shared_pct := sprintf('%0.2f%%',shared/tested*100)]
merged[,glucose_cbn := sprintf('%s (%s)',glucose,glucose_pct)]
merged[,galactose_cbn := sprintf('%s (%s)',galactose,galactose_pct)]
merged[,shared_cbn := sprintf('%s (%s)',shared,shared_pct)]
out_fn = sprintf('%s/table1.txt',out_dir)
fwrite(merged,out_fn,sep='\t')
