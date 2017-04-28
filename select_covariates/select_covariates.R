library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)

# Variables:
cov_dir='../processed_data/select_covariates/select_covariates/covariates/'
log_dir='../logs/select_covariates/'
out_dir='../processed_data/select_covariates/select_covariates/chr22/'

if (!dir.exists(cov_dir)) {dir.create(cov_dir,recursive=TRUE)}
if (!dir.exists(log_dir)) {dir.create(log_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Glucose:
cov_fn='../processed_data/select_covariates/merge_covariates/glu_sex_geno_sva.tsv'
cov=fread(cov_fn)
setorder(cov,sample)
cov[,sample:=NULL]

for (i in 1:ncol(cov)){
	x=cov[,1:i,with=F]

	xtxt=sprintf('%s/glu_cov.%s.txt',cov_dir,i)
	fwrite(x,xtxt,sep='\t',col.names=F,row.names=F)

	xbin=gsub('txt','bin',xtxt)
	fxbin=file(xbin,"wb")
	writeBin(as.double(c(as.matrix(x))), fxbin)

	sub_dir=sprintf("%s/cov%s",out_dir,i)
	if (!dir.exists(sub_dir)) {dir.create(sub_dir,recursive=TRUE)}
	tmp=foreach(j=1:50)%dopar%{
		command=sprintf("bash rasqual/rasqual.sh ../processed_data/rasqual/input/rasqual.input.chr22.filt.txt %s ../processed_data/rasqual/expression/glucose.expression.bin ../processed_data/rasqual/expression/glucose.size_factors_gc.bin ../data/genotype/asvcf/glucose_sid/rpe.imputed.chr22.all_filters.vcf.new.gz joint %s 2> %s/rasqual.chr22.%s.log",j,sub_dir,xbin,log_dir,j)
		system(command)
	}
}


# Galactose: 
# cov_fn='../processed_data/select_covariates/merge_covariates/gal_sex_geno_sva.tsv'
# cov=fread(cov_fn)
# setorder(cov,sample)
# cov[,sample:=NULL]

# for (i in 1:ncol(cov)){
# 	x=cov[,1:i,with=F]

# 	xtxt=sprintf('%s/gal_cov.%s.txt',cov_dir,i)
# 	fwrite(x,xtxt,sep='\t',col.names=F,row.names=F)

# 	xbin=gsub('txt','bin',xtxt)
# 	fxbin=file(xbin,"wb")
# 	writeBin(as.double(c(as.matrix(x))), fxbin)

# 	tmp=foreach(j=1:50)%dopar%{
# 		command=sprintf("bash rasqual/rasqual.sh ../processed_data/rasqual/input/rasqual.input.chr22.filt.txt %s ../processed_data/rasqual/expression/glucose.expression.bin ../processed_data/rasqual/expression/glucose.size_factors_gc.bin ../data/genotype/asvcf/glucose_sid/rpe.imputed.chr22.all_filters.vcf.new.gz joint %s %s 2> %s/rasqual.chr22.%s.log",j,out_dir,xbin,log_dir,j)
# 		system(command)
# 	}

# }
