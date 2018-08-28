# test_association.R tests association between SNP and metabolic phenotypes. 
library(data.table)
library(nlme)
library(doMC)
library(foreach)
registerDoMC(15)
library(gap)

phenotype_fn = '../processed_data/phenotype_association/format_phenotype/phenotypes.txt'
dosage_fn = '../processed_data/phenotype_association/format_genotype//dosage.txt'
covariate_fn = '../processed_data/select_covariates/merge_covariates/glu_sex_geno_sva.tsv'
fig_dir = '../figures/phenotype_association/test_association/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
out_dir = '../processed_data/phenotype_association/test_association/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

phenotype = fread(phenotype_fn)
dosage = fread(dosage_fn)
covariate = fread(covariate_fn)[,1:5]

phenotype_list = unique(phenotype$variable)
rsid_list = unique(dosage$rsid)

phenotype = phenotype[sample!='021011']

result_all = data.table()
for (pheno in phenotype_list){

	print(pheno)
	x = phenotype[variable == pheno,list(sample, value)]

	result = foreach(rs = rsid_list,.combine='rbind')%dopar%{

		print(rs)
		y = dosage[rsid == rs, list(sample, dosage)]

		tryCatch(expr = {
			data = merge(x,y,by='sample')
			data = merge(data,covariate,by='sample')
			data[, sample := as.factor(sample)]

			fit = lme(value ~ 1 + dosage + gender + PC1 + PC2 + PC3,
				random = ~ 1 | sample,
				data = data, method = 'REML')

			res = data.table(t(summary(fit)$tTable[2,]))
			res$pheno = pheno
			res$snp = rs
			return(res)
		},
			error = function(e){return(data.table())}
		)

	}

	result_all = rbind(result_all,result)

}

out_fn = sprintf('%s/all_association.txt', out_dir)
fwrite(result_all, out_fn, sep='\t')
fig_fn = sprintf('%s/qqplot.pdf',fig_dir)
pdf(fig_fn)
for (phe in phenotype_list){
	hist(result_all[pheno == phe,`p-value`],xlab='p-value',main = phe)
	qqunif(result_all[pheno == phe,`p-value`])
}
dev.off()


# Test glucose_bound_os and galactose_bound_os after
# removing the sample 020311. 

phenotype = phenotype[sample!='020311' & variable %in% c('glucose_bound_os','galactose_bound_os'),]
phenotype_list = unique(phenotype$variable)

result_all = data.table()
for (pheno in phenotype_list){

	print(pheno)
	x = phenotype[variable == pheno,list(sample, value)]

	result = foreach(rs = rsid_list,.combine='rbind')%dopar%{

		print(rs)
		y = dosage[rsid == rs, list(sample, dosage)]

		tryCatch(expr = {
			data = merge(x,y,by='sample')
			data = merge(data,covariate,by='sample')
			data[, sample := as.factor(sample)]

			fit = lme(value ~ 1 + dosage + gender + PC1 + PC2 + PC3,
				random = ~ 1 | sample,
				data = data, method = 'REML')

			res = data.table(t(summary(fit)$tTable[2,]))
			res$pheno = pheno
			res$snp = rs
			return(res)
		},
			error = function(e){return(data.table())}
		)

	}

	result_all = rbind(result_all,result)

}

out_fn = sprintf('%s/bound_os_no_020311_association.txt', out_dir)
fwrite(result_all, out_fn, sep='\t')
fig_fn = sprintf('%s/bound_os_no_020311_qqplot.pdf',fig_dir)
pdf(fig_fn)
for (phe in phenotype_list){
	hist(result_all[pheno == phe,`p-value`],xlab='p-value',main = phe)
	qqunif(result_all[pheno == phe,`p-value`])
}
dev.off()