Metabolic phenotype analysis find SNPs correlated with metabolic phenotypes. 

####################
# format phenotype #
####################


# Main
declare number of worksheets
declare output directory

for i in number of worksheets: 
	
	data = read excel sheets i

end for 

write data to output


# Functions
declare excel_fn = excel file name
declare phenotypes = column headers

method: read_excel_sheet(excel_fn, sheetIndex)
	
	read excel sheet 

	remove first two rows (headers)

	change the column name to be the phenotypes 

	melt the data.frame into 2 columns: {variable, value}

	remove NA rows

	add id column

	return melt_data

end method

###################
# format genotype #
###################
# Main: 

read eGenes for both conditions (we will test all lead SNPs, regardless of condition)

read lead SNPs from the glucose condition, breaking ties by selecting the first SNP. 

for snp in lead SNPs:

	data = read VCF at the snp location

end for 

convert sample id to RNA id

write data to output 

declare vcf_dir 
method: read VCF(vcf_dir, chr, pos):

	vcf_fn = find VCF file from vcf_dir

	bcftools to subset VCF 

	melt data 

	add SNP identifier

end method


###################
# Model selection #
###################

1. test full model with no random effect
2. test full model with random intercept 
3. test full model with random intercept and slope 

use AIC or BIC to select a model 

####################
# test association #
####################

declare phenotypes_fn
declare dosage_fn
declare covariate_fn 

read phenotype
read dosage 
covariate = read covariate


for pheno in phenotypes:

	for snp in SNP:

		data = merge pheno and snp

		data = merge data and covariante

		test association between snp and pheno

	end for 

end for

declare covariate_fn
method: read covarite(covariate_fn):

	read_covariate

end method

method: test association(data):

	test association between snp and pheno using mixed effects model while controlling for age, sex, and ethnicity.

end method


#################
# plot top hits #
#################
declare all_hits_fn 
declare pheontype_fn 
declare dosage_fn

read all_hits 
read phenotype 
read dosage 

convert dosage to genotype

calculate all_hits FDR 
subset to top_hits

for hit in top_hits:
	
	pheno = subset phenotype based on hit 

	geno = subset genotype based on hit 

	merged = merge pheno and geno

	plot merged

end for 


method: plot merged(merged):

	convert dosage to integer genotypes 

	plot genotype versus phenotype 

end method 