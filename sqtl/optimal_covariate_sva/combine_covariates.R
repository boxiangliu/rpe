library(data.table)
library(stringr)

glucose_sva_fn = '../processed_data/hidden_covariates/sva/v2/SVAFact_Spl_glucose_scaled_ZeroLeq12.txt'
galactose_sva_fn = '../processed_data/hidden_covariates/sva/v2/SVAFact_Spl_galactose_scaled_ZeroLeq12.txt'
gender_fn = '../processed_data/sex/sex/gender.tsv'
genotype_fn = '../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv'
out_dir = '../processed_data/optimal_covariate_sva/combine_covariates/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_sva = function(sva_fn){
	sva = fread(sva_fn)
	sva[, rn:=str_replace(rn,'_glucose','')]
	sva[, rn:=str_replace(rn,'_galactose','')]
	return(sva)
}

read_gender = function(gender_fn){
	gender = fread(gender_fn)
	gender[,gender:=ifelse(gender=='M',0,1)]
	return(gender)
}

read_genotype = function(genotype_fn){
	genotype = fread(genotype_fn)
	genotype = genotype[,list(sample,PC1,PC2,PC3)]
	return(genotype)
}

glucose_sva = read_sva(glucose_sva_fn)
galactose_sva = read_sva(galactose_sva_fn)
gender = read_gender(gender_fn)
genotype = read_genotype(genotype_fn)

merged = merge(gender,genotype,by='sample')

glucose_merged = merge(merged,glucose_sva,by.x='sample',by.y='rn')
glucose_merged = as.data.frame(t(glucose_merged))

fwrite(glucose_merged,sprintf('%s/glucose_covariate.txt',out_dir),col.names=FALSE,row.names=TRUE,sep='\t')

galactose_merged = merge(merged,galactose_sva,by.x='sample',by.y='rn')
galactose_merged = as.data.frame(t(galactose_merged))

fwrite(galactose_merged,sprintf('%s/galactose_covariate.txt',out_dir),col.names=FALSE,row.names=TRUE,sep='\t')
