library(data.table)
library(dplyr)
library(dtplyr)
library(stringr)

# Combine GTEx and RPE RPKM, 
# Perform log(x+2) transformation:

gtex_rpkm_dir='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_rpkm/'
rpe_rpkm_fn='../data/rnaseq/rpkm/rpe.rpkm'
tmp_dir='../processed_data/mds/tmp/'
out_dir='../processed_data/mds/preprocess/'
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

input_list=list.files(gtex_rpkm_dir,full.names=TRUE)
input_list=c(input_list,sprintf('%s/%s.rpkm.txt',tmp_dir,c('glucose','galactose')))


# read and combine rpkm files:
combined_rpkm=data.frame()
col_data=data.frame()
for (input in input_list){
	message('INFO - reading ',input)
	tissue=input %>% basename() %>% str_replace('.rpkm.txt','')

	# read a rpkm file: 
	rpkm=fread(input,header=T)
	stopifnot(length(unique(rpkm$Name))==length(rpkm$Name)) # sanity check
	setnames(rpkm,'Gene','Name')

	# if this rpkm file is the first being read: 
	if (ncol(combined_rpkm)==0){
		row_data=rpkm$Name
		rpkm[,Name:=NULL] 
		combined_rpkm=rpkm
		col_data=data.frame(sample=colnames(rpkm),tissue=tissue)
	} 
	# if this rpkm file is not the first being read: 
	else {
		# make the row orders of this rpkm dataframe consistent with previous rpkm dataframes: 
		if (!all(rpkm$Name==row_data)) {
			stopifnot(setequal(rpkm$Name,row_data))
			idx=match(row_data,rpkm$Name)
			rpkm=rpkm[idx,]
		}
		# merge the rpkm dataframe with previous rpkm dataframes:
		stopifnot(rpkm$Name==row_data)
		rpkm[,Name:=NULL]
		combined_rpkm=cbind(combined_rpkm,rpkm)
		col_data=rbind(col_data, data.frame(sample=colnames(rpkm),tissue=tissue))
	}
}

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
combined_rpkm$Name=row_data
setcolorder(combined_rpkm,c(ncol(combined_rpkm),seq(ncol(combined_rpkm)-1)))


# write output: 
out_rpkm=paste0(out_dir,'/combined.rpkm')
out_col_data=paste0(out_dir,'/combined.col')
message('writing output...')
fwrite(combined_rpkm,out_rpkm,sep="\t")
fwrite(col_data,out_col_data,sep="\t")
