library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)

# Combine GTEx and RPE RPKM, 
# Perform log(x+2) transformation:

gtex_rpkm_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/rna-seq/GTEx_Analysis_v7_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz'
gtex_sample_info_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/sample_annotations/GTEx_Analysis_2016-01-15_v7_SampleAttributesDS.txt'
rpe_rpkm_fn='../data/rnaseq/rpkm/rpe.rnaseqc_rpkm'
tmp_dir='../processed_data/mds/tmp/'
out_dir='../processed_data/mds/preprocess.GTExV7/'
if (!dir.exists(tmp_dir)) dir.create(tmp_dir,recursive=TRUE)
if (!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)
rpkm_threshold=0.1
num_indv_threshold=10

# split RPE into glucose and galactose conditions: 
rpe=fread(rpe_rpkm_fn)


glucose=rpe%>%select(Name,contains('glucose'))
fwrite(glucose,sprintf('%s/glucose.rpkm.txt',tmp_dir),sep='\t')

galactose=rpe%>%select(Name,contains('galactose'))
fwrite(galactose,sprintf('%s/galactose.rpkm.txt',tmp_dir),sep='\t')

gtex=fread(paste('gunzip -c', gtex_rpkm_fn))
gtex_sample_info=fread(gtex_sample_info_fn,select=c(1,14))
setnames(gtex_sample_info,c('sample','tissue'))


# read and combine rpkm files:
combined_rpkm=merge(gtex,rpe,by=c('Name','Description'))
row_data=combined_rpkm[,list(Name,Description)]

setDF(combined_rpkm)
rownames(combined_rpkm)=combined_rpkm$Name
combined_rpkm$Name=NULL
combined_rpkm$Description=NULL

temp=data.table(sample=colnames(rpe)[3:ncol(rpe)])
temp[,tissue:=ifelse(str_detect(sample,'lucose'),'RPE (glu)', 'RPE (gal)')]


col_data=rbind(gtex_sample_info,temp)
idx=match(colnames(combined_rpkm),col_data$sample)
col_data=col_data[idx,]

stopifnot(col_data$sample==colnames(combined_rpkm))


# filter for genes with rpkm>0.1 in > 10 individuals:
message('filtering...')
pass=combined_rpkm>rpkm_threshold
keep=rowSums(pass)>num_indv_threshold
combined_rpkm=combined_rpkm[keep,]
row_data=row_data[keep]


# log2(x+2) transform:
message('log transforming...')
combined_rpkm_log2xp1=log2(combined_rpkm+1)
combined_rpkm_log2xp2=log2(combined_rpkm+2)


# add gene id: 
combined_rpkm$Name=row_data$Name
setcolorder(combined_rpkm,c(ncol(combined_rpkm),seq(ncol(combined_rpkm)-1)))
combined_rpkm_log2xp1$Name=row_data$Name
setcolorder(combined_rpkm_log2xp1,c(ncol(combined_rpkm_log2xp1),seq(ncol(combined_rpkm_log2xp1)-1)))
combined_rpkm_log2xp2$Name=row_data$Name
setcolorder(combined_rpkm_log2xp2,c(ncol(combined_rpkm_log2xp2),seq(ncol(combined_rpkm_log2xp2)-1)))


# write output: 
message('writing output...')
out_rpkm=paste0(out_dir,'/combined.rpkm')
out_col_data=paste0(out_dir,'/combined.col')
fwrite(combined_rpkm,out_rpkm,sep="\t")
fwrite(col_data,out_col_data,sep="\t")

out_rpkm=paste0(out_dir,'/combined_logxp1.rpkm')
fwrite(combined_rpkm_log2xp1,out_rpkm,sep="\t")

out_rpkm=paste0(out_dir,'/combined_logxp2.rpkm')
fwrite(combined_rpkm_log2xp2,out_rpkm,sep="\t")

