library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)
library(gap)
library(stringr)


annotation_fn = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6_v6p_annotations/gencode.v19.genes.v6p.patched_contigs.bed'
fig_dir = '../figures/gwas_QTL_overlap/qqplot/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

amd_gwas_fn = '../data/gwas/Fritsche_2015_AdvancedAMD.txt'
myopia_gwas_fn = '../data/gwas/23andme_myopia.prepared.txt.gz'

MT_treeQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
glu_rasqual_dir = '../processed_data/rasqual/output/glucose/joint/'
gal_rasqual_dir = '../processed_data/rasqual/output/galactose/joint/'

glu_sqtl_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
gal_sqtl_fn = '../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz'

read_treeQTL_eGenes = function(fn,annotation_fn){
	x = fread(fn)[,list(gene_id = gene)]
	annotation = fread(annotation_fn)[,c(1,13)][,list(
		chr = paste0('chr',V1), gene_id = V13)]
	x = merge(x,annotation,by='gene_id')
	return(x)
}

read_MT_treeQTL_eGenes = function(fn,annotation_fn){
	x = fread(fn)
	x = x[glucose == 1 & galactose == 1,list(gene_id = gene)]
	annotation = fread(annotation_fn)[,c(1,13)][,list(
		chr = paste0('chr',V1), gene_id = V13)]
	x = merge(x,annotation,by='gene_id')
	return(x)
}

parse_SNP = function(snp){
	parsed_snp = str_split_fixed(snp,'_',5)
	chr = paste0('chr',parsed_snp[,1])
	pos = as.integer(parsed_snp[,2])
	return(list(chr,pos))
}

read_sQTL = function(fn,fdr = 0.05){
	if (grepl('.gz',fn)){
		x = fread(sprintf('gunzip -c %s',fn))
	} else {
		x = fread(fn)
	}
	y = x[,c(1,6,16)]
	setnames(y,c('cluster','snp','pval'))
	y[,padj := p.adjust(pval,method="fdr")]
	z = y[padj <= fdr]
	z[,c('chr','pos'):=parse_SNP(snp)]
	return(z)
}

eGenes_to_eSNPs = function(eGenes,rasqual_dir){
	eSNPs = foreach(i = 1:nrow(eGenes), .combine = 'rbind')%dopar%{
		g = eGenes[i,gene_id]
		chr = eGenes[i,chr]
		dir = paste(rasqual_dir,chr,sep='/')
		fn = list.files(dir, pattern = g, full.names = TRUE)
		stopifnot(length(fn) == 1)
		x = fread(fn)[,c(1,2,3,4,11)][,list(gene = V1, 
			snp = V2, chr = V3, pos = V4, chisq = V11)]
		snp = x[which.max(chisq),]
		snp[,gene := str_split_fixed(gene,'_',2)[,1]]
		return(snp)
	}
	return(eSNPs)
}

eGenes_to_GTEx_eSNPs = function(top=nrow(eGenes),exclude,GTEx_top_association_fn){
	command = paste('zcat',GTEx_top_association_fn)
	x = fread(command)[,list(gene=gene_id,snp=variant_id,chr,pos,pval_beta)]
	x = x[!(gene %in% exclude),]
	x[,chr:=paste0('chr',chr)]
	x[,rank:=rank(pval_beta)]
	eSNPs = x[rank<=top,]
	return(eSNPs)
}

read_GWAS = function(fn,rsid_col,chr_col,pos_col,pval_col){
	if (grepl('.gz',fn)){
		gwas = fread(sprintf('gunzip -c %s',fn))
	} else {
		gwas = fread(fn)
	}
	setnames(gwas,c(rsid_col,chr_col,pos_col,pval_col),
		c('rsid','chr','pos','pval'))
	if (!grepl('chr',gwas$chr[1])){
		gwas[,chr := paste0('chr',chr)]
	}
	gwas = gwas[,list(rsid,chr,pos,pval)]
	return(gwas)
}

select_bg_SNP = function(){
	return()
}

subset_GWAS = function(gwas, subset){
	stopifnot(c('chr','pos') %in% colnames(gwas))
	stopifnot(c('chr','pos') %in% colnames(subset))

	subset = merge(gwas,subset[,list(chr,pos)],by=c('chr','pos'))
	return(subset)
}

make_qqplot = function(pval_list,pval_threshold = 1e-16, n = 1e5){
	p = pval_list[[1]]
	if (length(p) > n){
		p = sample(p, n)
	}
	p[p < pval_threshold] = pval_threshold
	qqunif(p, col = 1, pch = 16)

	for (i in 2:length(pval_list)){
		p = pval_list[[i]]
		if (length(p) > n){
			p = sample(p, n)
		}
		p[p < pval_threshold] = pval_threshold
		points = qqunif(p, plot.it = FALSE)
		points(points, col = i, pch = 16)
	}

	legend('topleft',pch = 16, col = 1:length(p), 
		legend = names(pval_list))
}

add_GTEx_eSNP_to_pval_list = function(pval_list,gwas){
	LCL_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Cells_EBV-transformed_lymphocytes.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,LCL_eSNPs)
	pval_list[['LCL']] = subset$pval

	blood_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Whole_Blood.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,blood_eSNPs)
	pval_list[['blood']] = subset$pval

	skin_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Skin_Sun_Exposed_Lower_leg.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,skin_eSNPs)
	pval_list[['skin']] = subset$pval

	return(pval_list)
}

amd_gwas = read_GWAS(amd_gwas_fn, rsid_col = 'Marker', 
	chr_col = 'Chrom', pos_col = 'Pos', pval_col = 'GC.Pvalue')


pval_list = list(GWAS = amd_gwas$pval)
eQTL_fn_list = c(`Glucose eQTL`=glu_rasqual_dir,`Galactose eQTL`=gal_rasqual_dir)
for (i in seq_along(eQTL_fn_list)){
	rasqual_dir = eQTL_fn_list[i]
	name = names(eQTL_fn_list)[i]
	eGenes = read_MT_treeQTL_eGenes(MT_treeQTL_fn,annotation_fn)
	eSNPs = eGenes_to_eSNPs(eGenes,rasqual_dir)
	subset = subset_GWAS(amd_gwas,eSNPs)
	pval_list[[name]] = subset$pval
}

pval_list = add_GTEx_eSNP_to_pval_list(pval_list,amd_gwas)

pdf(sprintf('%s/%s_GWAS_eQTL_qqplot_GTEx_eQTL.pdf',fig_dir,'AMD'))
make_qqplot(pval_list)
dev.off()

myopia_gwas = read_GWAS(myopia_gwas_fn, rsid_col = 'rsid',
	chr_col = 'chr', pos_col = 'snp_pos', pval_col = 'pvalue')

pval_list = list(GWAS = myopia_gwas$pval)
eQTL_fn_list = c(`Glucose eQTL`=glu_rasqual_dir,`Galactose eQTL`=gal_rasqual_dir)
for (i in seq_along(eQTL_fn_list)){
	rasqual_dir = eQTL_fn_list[i]
	name = names(eQTL_fn_list)[i]
	eGenes = read_MT_treeQTL_eGenes(MT_treeQTL_fn,annotation_fn)
	eSNPs = eGenes_to_eSNPs(eGenes,rasqual_dir)
	subset = subset_GWAS(myopia_gwas,eSNPs)
	pval_list[[name]] = subset$pval
}

pval_list = add_GTEx_eSNP_to_pval_list(pval_list,myopia_gwas)
pdf(sprintf('%s/%s_GWAS_eQTL_qqplot_GTEx_eQTL.pdf',fig_dir,'Myopia'))
make_qqplot(pval_list)
dev.off()


