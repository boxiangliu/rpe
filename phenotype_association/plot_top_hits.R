# Plot top hits makes scatter plot for the top geno-pheno associations

library(data.table)
library(cowplot)

all_hits_fn = '../processed_data/phenotype_association/test_association/all_association.txt'
phenotype_fn = '../processed_data/phenotype_association/format_phenotype/phenotypes.txt'
dosage_fn = '../processed_data/phenotype_association/format_genotype//dosage.txt'
fig_dir = '../figures/phenotype_association/plot_top_hits/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}


plot_data = function(data){
	ggplot(data, aes(x = as.factor(genotype), y = value)) + 
		geom_boxplot() + 
		geom_jitter(aes(color = sample)) + 
		xlab('Genotype') + 
		ylab('Measurement')
}

# read data
all_hits = fread(all_hits_fn)
phenotype = fread(phenotype_fn)
phenotype = phenotype[sample!='021011']
dosage = fread(dosage_fn)

# convert dosage to genotype
dosage$genotype = round(dosage$dosage)

# calculate all_hits FDR
all_hits[,fdr:=p.adjust(`p-value`,method = 'fdr')]
sum(all_hits[,fdr] < 0.05) # 0

# plot: 
setorder(all_hits,`p-value`)
fig_fn = sprintf('%s/top_hits.pdf',fig_dir)
pdf(fig_fn)
for (i in 1:5){
	r = all_hits[i,snp]
	ph = all_hits[i,pheno]

	pheno = phenotype[variable == ph]
	snp = dosage[rsid == r]

	data = merge(pheno, snp, by = 'sample')
	p = plot_data(data)
	p = p + ylab(ph) + xlab(r)
	print(p)

}
dev.off()

# plot again after removing 020311 (outlier):
fig_fn = sprintf('%s/top_hits_no_020311.pdf',fig_dir)
pdf(fig_fn)
for (i in 1:5){

	r = all_hits[i,snp]
	ph = all_hits[i,pheno]

	pheno = phenotype[variable == ph]
	snp = dosage[rsid == r]

	data = merge(pheno, snp, by = 'sample')
	data_2 = data[sample!='020311']
	p_2 = plot_data(data_2)
	p_2 = p_2 + ylab(ph) + xlab(r)
	print(p_2)

}
dev.off()



