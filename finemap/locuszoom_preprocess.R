library(data.table)
library(stringr)

gwas_fn='../data/gwas/Fritsche_2015_AdvancedAMD.txt'
eqtl_fn='../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135437.5_RDH5.txt'
out_dir='../processed_data/finemap/locuszoom_preprocess/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

gwas=fread(gwas_fn)
eqtl=fread(eqtl_fn,select=c(1,3,4,5,6,11),col.names=c('gene','chr','pos','ref','alt','chisq'))
eqtl[,pval:=pchisq(chisq,df=1,lower.tail=FALSE)]

eqtl[,chr:=str_replace(chr,'chr','')]

eqtl=merge(eqtl,gwas[,list(chr=Chrom,pos=Pos,Marker)],by=c('chr','pos'))
fwrite(eqtl,sprintf('%s/RHD5.txt',out_dir),sep='\t')