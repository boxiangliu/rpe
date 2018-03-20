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

glu_treeQTL_fn = '../processed_data/rasqual/output/glucose/treeQTL/eGenes.txt'
gal_treeQTL_fn = '../processed_data/rasqual/output/galactose/treeQTL/eGenes.txt'
glu_rasqual_dir = '../processed_data/rasqual/output/glucose/joint/'
gal_rasqual_dir = '../processed_data/rasqual/output/galactose/joint/'

glu_sqtl_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
gal_sqtl_fn = '../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz'

read_treeQTL_eGenes = function(fn,annotation_fn){
	x = fread(fn)[,list(gene_id = family)]
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
		return(snp)
	}
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



amd_gwas = read_GWAS(amd_gwas_fn, rsid_col = 'Marker', 
	chr_col = 'Chrom', pos_col = 'Pos', pval_col = 'GC.Pvalue')


pval_list = list(GWAS = amd_gwas$pval)
eQTL_fn_list = list(`Glucose eQTL`=c(glu_treeQTL_fn,glu_rasqual_dir),
	`Galactose eQTL`=c(gal_treeQTL_fn,gal_rasqual_dir))
for (i in seq_along(eQTL_fn_list)){
	fn = eQTL_fn_list[[i]][1]
	rasqual_dir = eQTL_fn_list[[i]][2]
	name = names(eQTL_fn_list)[i]
	eGenes = read_treeQTL_eGenes(fn,annotation_fn)
	eSNPs = eGenes_to_eSNPs(eGenes,rasqual_dir)
	subset = subset_GWAS(amd_gwas,eSNPs)
	pval_list[[name]] = subset$pval
}

sQTL_fn_list = c(`Glucose sQTL` = glu_sqtl_fn,
	`Galactose sQTL` = gal_sqtl_fn)
for (i in seq_along(sQTL_fn_list)){
	fn = sQTL_fn_list[i]
	name = names(sQTL_fn_list)[i]
	sSNPs = read_sQTL(fn)
	subset = subset_GWAS(amd_gwas,sSNPs)
	pval_list[[name]] = subset$pval
}


pdf(sprintf('%s/%s_GWAS_eQTL_qqplot.pdf',fig_dir,'AMD'))
make_qqplot(pval_list)
dev.off()

















