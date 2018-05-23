library(data.table)
amd_gwas_fn = '../data/gwas/Fritsche_2015_AdvancedAMD.txt'
out_fn = '../data/gwas/Fritsche_2015_AdvancedAMD.v2.txt'
amd_gwas = fread(amd_gwas_fn)
amd_gwas[,effect_allele := ifelse(Overall=='+',Allele1,Allele2)]
amd_gwas[,non_effect_allele := ifelse(Overall=='-',Allele1,Allele2)]
amd_gwas[,c('Allele1','Allele2','Overall'):=NULL]

fwrite(amd_gwas,out_fn,sep='\t')