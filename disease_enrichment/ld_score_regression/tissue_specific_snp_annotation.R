library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)

# Variables:
tissue_specific_gene_dir='../processed_data/rpe_specific_genes.GTExV7/'
kept_tissue_fn='../processed_data/rpe_specific_genes.GTExV7/tissue_kept.txt'
gene_annotation_fn='../data/reference/gencode.v19.annotation.collapsed_annotation.gtf'
plink_dir='/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/'
out_dir='../processed_data/disease_enrichment/ld_score_regression/tissue_specific_snp_annotation/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Functions:
read_tissue_specific_gene=function(tissue_specific_gene_fn_list,kept_tissue=NULL){
	tissue_specific_gene=foreach(i=seq_along(tissue_specific_gene_fn_list),.combine='rbind')%dopar%{
		tissue_specific_gene_fn=tissue_specific_gene_fn_list[i]
		tissue = str_split_fixed(basename(tissue_specific_gene_fn),'\\.',3)[,2]
		tissue = str_replace_all(tissue,'_',' ')

		if (tissue %in% kept_tissue | is.null(kept_tissue)){
			data.table(
				fread(tissue_specific_gene_fn,select=1,col.names='GENE_ID'),
				TISSUE=tissue)
		} else {
			data.table()
		}

	}
	return(tissue_specific_gene)
}

read_top_tissue_specific_gene=function(tissue_specific_gene_fn_list,top,kept_tissue=NULL){
	tissue_specific_gene=foreach(i=seq_along(tissue_specific_gene_fn_list),.combine='rbind')%dopar%{
		tissue_specific_gene_fn=tissue_specific_gene_fn_list[i]
		tissue = str_split_fixed(basename(tissue_specific_gene_fn),'\\.',3)[,2]
		tissue = str_replace_all(tissue,'_',' ')

		if (tissue %in% kept_tissue | is.null(kept_tissue)){
			x = data.table(
				fread(tissue_specific_gene_fn,select=c(1,5,6,11)),
				TISSUE=tissue)
			setnames(x,'gene_id','GENE_ID')
			x = x[type=='protein_coding'&chr%in%paste0('chr',1:22)&!is.na(zscore),]
			setorder(x,-zscore)
			x = x[1:top,list(GENE_ID,TISSUE)]
		} else {
			x = data.table()
		}
		return(x)
	}
	return(tissue_specific_gene)
}

read_gene_annotation = function(gene_annotation_fn){
	gene_annotation=fread(gene_annotation_fn,
		select=c(1,3:5,9),
		col.names=c('CHR','TYPE','START','END','ANNOT'))
	gene_annotation[,GENE_ID:=str_extract(ANNOT,'(?<=gene_id \\")(ENSG[0-9\\.]+)(?=\\";)')]
	gene_annotation[,ANNOT:=NULL]
}

merge_exon_and_tissue_specific_gene=function(gene_annotation,tissue_specific_gene){
	tissue_specific_exon_interval=merge(gene_annotation[TYPE=='exon',],tissue_specific_gene,by='GENE_ID',allow.cartesian=TRUE)
	setorder(tissue_specific_exon_interval,TISSUE,CHR,START)
	stopifnot(tissue_specific_exon_interval[,all(END>=START)])
	tissue_specific_exon_interval[,c('START','END'):=list(START-1000,END+1000)]
	return(tissue_specific_exon_interval)
}

read_bim = function(plink_dir){
	bim_fn_list=list.files(plink_dir,'bim',full.names=TRUE)
	bim=foreach(i=seq_along(bim_fn_list),.combine='rbind')%dopar%{
		bim_fn=bim_fn_list[i]
		message(bim_fn)
		fread(bim_fn,col.names=c('CHR','SNP','CM','BP','A1','A2'))
	}
	bim[,c('CHR','START','END'):=list(as.character(CHR),BP,BP)]
	setkey(bim,CHR,START,END)
	return(bim)
}

merge_tissue_specific_exon_and_bim=function(tissue_specific_exon_interval_wide,bim,tissue_list){
	setkey(tissue_specific_exon_interval_wide,CHR,START,END)
	temp=foverlaps(bim,tissue_specific_exon_interval_wide)
	annot=unique(temp[,c('CHR','BP','SNP','CM',tissue_list),with=F])
	
	setDF(annot)
	annot[is.na(annot)]=0
	
	setDT(annot)
	dup=which(duplicated(annot[,list(CHR,BP,SNP,CM)]))

	for (tissue in tissue_list){
		setnames(annot,tissue,'temp')
		annot[c(dup,dup-1),temp:=max(temp,na.rm=TRUE),by=c('CHR','BP','SNP','CM')]
		setnames(annot,'temp',tissue)
	}

	annot=unique(annot)
	stopifnot(annot$SNP==bim$SNP)

	return(annot)
}

output_merged_annotations=function(merged_annot,out_dir){
	if(!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
	foreach(chr=seq(22))%dopar%{
		message(chr)
		out_fn=gzfile(sprintf('%s/merged.%s.annot.gz',out_dir,chr))
		write.table(merged_annot[CHR==chr],out_fn,sep='\t',quote=FALSE,col.names=T,row.names=F,na='0')
	}
}

# Read and merge PLINK bim files:
bim = read_bim(plink_dir)

# Read gene annotation:
gene_annotation = read_gene_annotation(gene_annotation_fn)


# Read kept tissues:
kept_tissue=NULL

# Read tissue-specific genes:
message('INFO - Reading tissue-specific gene')
tissue_specific_gene_fn_list = list.files(tissue_specific_gene_dir,pattern='specific_genes',full.names=TRUE)
tissue_specific_gene = read_tissue_specific_gene(tissue_specific_gene_fn_list)


# Merge exon intervals and tissue specific genes:
message('INFO - Merging exon interval and tissue specific genes')
tissue_specific_exon_interval=merge_exon_and_tissue_specific_gene(gene_annotation,tissue_specific_gene)
tissue_specific_exon_interval_wide=dcast(tissue_specific_exon_interval,CHR+START+END~TISSUE,value.var='GENE_ID',fun.aggregate=length)
tissue_specific_exon_interval_wide[,CHR:=str_replace(CHR,'chr','')]


# Merge tissue-specific exon interval and bim file:
message('INFO - Merging tissue-specific exon interval and bim file')
tissue_list=sort(unique(tissue_specific_exon_interval$TISSUE))
annot=merge_tissue_specific_exon_and_bim(tissue_specific_exon_interval_wide,bim,tissue_list)
setnames(annot,names(annot),str_replace_all(names(annot),' ','_'))

# Output merged annotation by chromosome:
message('INFO - Outputing merged annotation by chromosome')
out_dir_2 = sprintf('%s/4sd/',out_dir)
if (!dir.exists(out_dir_2)) {dir.create(out_dir_2,recursive=TRUE)}
output_merged_annotations(annot,out_dir_2)


# Do the same for top tissue-specific genes:
tissue_specific_gene_fn_list = list.files(tissue_specific_gene_dir,pattern='all_genes',full.names=TRUE)
foreach(N = c(200,500,1000,2000,4000)) %dopar% {
	tissue_specific_gene = read_top_tissue_specific_gene(tissue_specific_gene_fn_list,top=N)

	# Merge exon intervals and tissue specific genes:
	message('INFO - Merging exon interval and tissue specific genes')
	tissue_specific_exon_interval=merge_exon_and_tissue_specific_gene(gene_annotation,tissue_specific_gene)
	tissue_specific_exon_interval_wide=dcast(tissue_specific_exon_interval,CHR+START+END~TISSUE,value.var='GENE_ID',fun.aggregate=length)
	tissue_specific_exon_interval_wide[,CHR:=str_replace(CHR,'chr','')]


	# Merge tissue-specific exon interval and bim file:
	message('INFO - Merging tissue-specific exon interval and bim file')
	tissue_list=sort(unique(tissue_specific_exon_interval$TISSUE))
	annot=merge_tissue_specific_exon_and_bim(tissue_specific_exon_interval_wide,bim,tissue_list)
	setnames(annot,names(annot),str_replace_all(names(annot),' ','_'))

	# Output merged annotation by chromosome:
	message('INFO - Outputing merged annotation by chromosome')
	out_dir_2 = sprintf('%s/top%s/',out_dir,N)
	if (!dir.exists(out_dir_2)) {dir.create(out_dir_2,recursive=TRUE)}
	output_merged_annotations(annot,out_dir_2)
}