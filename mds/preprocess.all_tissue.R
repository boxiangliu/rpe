library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)

# Combine GTEx and RPE RPKM, 
# Perform log(x+2) transformation:

gtex_rpkm_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_All_Tissues_gene_rpkm_FOR_QC_ONLY.gct'
gtex_sample_info_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/sample_annotations/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt'
rpe_rpkm_fn='../data/rnaseq/rpkm/rpe.rpkm'
tmp_dir='../processed_data/mds/tmp/'
out_dir='../processed_data/mds/preprocess.all_tissue/'
dir.create(tmp_dir,recursive=TRUE);on.exit(unlink(tmp_dir,recursive=TRUE))
dir.create(out_dir,recursive=TRUE)
rpkm_threshold=0.1
num_indv_threshold=10

# split RPE into glucose and galactose conditions: 
rpe=fread(rpe_rpkm_fn)
setnames(rpe,'gene_id','Gene')

glucose=rpe%>%select(Gene,contains('glucose'))
fwrite(glucose,sprintf('%s/glucose.rpkm.txt',tmp_dir),sep='\t')

galactose=rpe%>%select(Gene,contains('galactose'))
fwrite(galactose,sprintf('%s/galactose.rpkm.txt',tmp_dir),sep='\t')

gtex=fread(gtex_rpkm_fn)
gtex_sample_info=fread(gtex_sample_info_fn,select=c(1,14))
setnames(gtex_sample_info,c('sample','tissue'))


# read and combine rpkm files:
combined_rpkm=merge(gtex,rpe,by.x='Name',by.y='Gene')

row_data=combined_rpkm[,list(Name,Description)]

setDF(combined_rpkm)
rownames(combined_rpkm)=combined_rpkm$Name
combined_rpkm$Name=NULL
combined_rpkm$Description=NULL

temp=data.table(sample=colnames(rpe)[2:ncol(rpe)])
temp[,tissue:=ifelse(str_detect(sample,'glucose'),'RPE (glu)', 'RPE (gal)')]

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
combined_rpkm=log2(combined_rpkm+2)


# add gene id: 
combined_rpkm$Name=row_data$Name
setcolorder(combined_rpkm,c(ncol(combined_rpkm),seq(ncol(combined_rpkm)-1)))


# write output: 
out_rpkm=paste0(out_dir,'/combined.rpkm')
out_col_data=paste0(out_dir,'/combined.col')
message('writing output...')
fwrite(combined_rpkm,out_rpkm,sep="\t")
fwrite(col_data,out_col_data,sep="\t")
