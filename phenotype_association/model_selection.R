# test_association.R tests association between SNP and metabolic phenotypes. 
library(data.table)
library(nlme)

phenotype_fn = '../processed_data/phenotype_association/format_phenotype/phenotypes.txt'
dosage_fn = '../processed_data/phenotype_association/format_genotype//dosage.txt'

phenotype = fread(phenotype_fn)
dosage = fread(dosage_fn)


phenotype_list = unique(phenotype$variable)
rsid_list = unique(dosage$rsid)

phenotype = phenotype[sample!='021011']
pheno = phenotype_list[1]
rs = rsid_list[1]

x = phenotype[variable == pheno,list(sample, value)]
y = dosage[rsid == rs, list(sample, dosage)]
data = merge(x,y,by='sample')
data[, sample := as.factor(sample)]

b1 = gls(value ~ 1 + dosage, data = data, method = 'REML')
b2 = lme(value ~ 1 + dosage,
	random = ~ 1 | sample,
	data = data, method = 'REML')
b3 = lme(value ~ 1 + dosage,
	random = ~ 1 + dosage | sample,
	data = data, method = 'REML')

AIC(b1,b2,b3) # b3 is best
BIC(b1,b2,b3) # b2 is best 
anova(b1,b2,b3) # b3 is best
0.5*(1-pchisq(95.65992,1)) # 0
0.5*((1-pchisq(5.00418,1)) + (1-pchisq(5.00418,1))) # 0.025

plot(b3$fitted[,'fixed'],b3$residuals[,'fixed'])