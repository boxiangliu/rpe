library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)
library(gap)
library(stringr)
source('utils/gtex_tissue_info.R')
library(cowplot)

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

tissue_color = c(tissue_color,GWAS='#000000')
tissue_abbreviation = c(tissue_abbreviation,GWAS='GWAS')

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

make_pval_list = function(gwas){
	pval_list = list(GWAS = gwas$pval)
	eQTL_fn_list = c(`RPE - glucose`=glu_rasqual_dir,`RPE - galactose`=gal_rasqual_dir)
	for (i in seq_along(eQTL_fn_list)){
		rasqual_dir = eQTL_fn_list[i]
		name = names(eQTL_fn_list)[i]
		eGenes = read_MT_treeQTL_eGenes(MT_treeQTL_fn,annotation_fn)
		eSNPs = eGenes_to_eSNPs(eGenes,rasqual_dir)
		subset = subset_GWAS(gwas,eSNPs)
		pval_list[[name]] = subset$pval
	}
	return(pval_list)
}

make_sqtl_pval_list = function(gwas){
	pval_list = list()
	sQTL_fn_list = c(`RPE - glucose` = glu_sqtl_fn,`RPE - galactose` = gal_sqtl_fn)
	for (i in seq_along(sQTL_fn_list)){
		fn = sQTL_fn_list[i]
		name = names(sQTL_fn_list)[i]
		sSNPs = read_sQTL(fn)
		subset = subset_GWAS(amd_gwas,sSNPs)
		pval_list[[name]] = subset$pval
	}
	return(pval_list)
}

subset_GWAS = function(gwas, subset){
	stopifnot(c('chr','pos') %in% colnames(gwas))
	stopifnot(c('chr','pos') %in% colnames(subset))

	subset = merge(gwas,subset[,list(chr,pos)],by=c('chr','pos'))
	return(subset)
}

add_GTEx_eSNP_to_pval_list = function(pval_list,gwas){
	LCL_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Cells_EBV-transformed_lymphocytes.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,LCL_eSNPs)
	pval_list[['Cells_EBV-transformed_lymphocytes']] = subset$pval

	blood_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Brain_Frontal_Cortex_BA9.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,blood_eSNPs)
	pval_list[['Brain_Frontal_Cortex_BA9']] = subset$pval

	skin_eSNPs = eGenes_to_GTEx_eSNPs(top=nrow(eGenes),eGenes$gene_id,'/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/top_associations/Skin_Sun_Exposed_Lower_leg.v7.egenes.txt.gz')
	subset = subset_GWAS(gwas,skin_eSNPs)
	pval_list[['Skin_Sun_Exposed_Lower_leg']] = subset$pval

	return(pval_list)
}

make_qqplot_data = function(pval_list,pval_threshold = 1e-16, n = 1e5){
	plot.new()
	x = foreach(i = 1:length(pval_list),.combine='rbind')%do%{
		p = pval_list[[i]]
		if (length(p) > n){
			p = sample(p, n)
		}
		p[p < pval_threshold] = pval_threshold
		x = data.frame(qqunif(p, plot.it = FALSE))
		x$tissue = names(pval_list)[i]
		return(x)
	}
	setDT(x)
	return(x)
}


plot_qqplot = function(qqplot_data){
	p = ggplot(qqplot_data,aes(x,y,color=tissue)) + 
		geom_point() + 
		geom_abline(intercept=0,slope=1,color='red') + 
		scale_color_manual(name='',values=tissue_color,labels = tissue_abbreviation) + 
		xlab(expression(-log[10]*'(expected p-value)')) + 
		ylab(expression(-log[10]*'(observed p-value)')) +
		theme(legend.position=c(0,1.1),legend.justification=c('left','top'))
	return(p)
}

plot_joint_qqplot = function(qqplot_data){
	p = ggplot(qqplot_data,aes(x,y,color=tissue,shape=type)) + 
		geom_point() + 
		geom_abline(intercept=0,slope=1,color='red') + 
		scale_color_manual(name='SNP set',values=tissue_color,labels = tissue_abbreviation) + 
		scale_shape_discrete(name='Phenotype') + 
		xlab(expression(-log[10]*'(expected p-value)')) + 
		ylab(expression(-log[10]*'(observed p-value)')) + 
		theme(legend.position=c(0,1),legend.justification=c('left','top'))
	return(p)
}
########
# eQTL #
########
# AMD:
amd_gwas = read_GWAS(amd_gwas_fn, rsid_col = 'Marker', 
	chr_col = 'Chrom', pos_col = 'Pos', pval_col = 'GC.Pvalue')

amd_pval_list = make_pval_list(amd_gwas)
amd_pval_list = add_GTEx_eSNP_to_pval_list(amd_pval_list,amd_gwas)
amd_qqplot_data = make_qqplot_data(amd_pval_list)
p1 = plot_qqplot(amd_qqplot_data)


# Myopia:
myopia_gwas = read_GWAS(myopia_gwas_fn, rsid_col = 'rsid',
	chr_col = 'chr', pos_col = 'snp_pos', pval_col = 'pvalue')

myopia_pval_list = make_pval_list(myopia_gwas)
myopia_pval_list = add_GTEx_eSNP_to_pval_list(myopia_pval_list,myopia_gwas)
myopia_qqplot_data = make_qqplot_data(myopia_pval_list)
p2 = plot_qqplot(myopia_qqplot_data)


# make plot:
p = plot_grid(p1+ggtitle('AMD'),p2+ggtitle('Myopia'),nrow=1,labels=c('A','B'))
fig_fn = sprintf('%s/RPE_GTEx_GWAS_qqplot.pdf',fig_dir)
save_plot(fig_fn,p,base_height=4,base_width=9)

fig_fn = sprintf('%s/RPE_GTEx_GWAS_qqplot.png',fig_dir)
save_plot(fig_fn,p,base_height=4,base_width=9)

########
# sQTL #
########
amd_sqtl_pval_list = make_sqtl_pval_list(amd_gwas)
amd_sqtl_qqplot_data = make_qqplot_data(amd_sqtl_pval_list)
amd_qqplot_data$type = 'eQTL'
amd_sqtl_qqplot_data$type = 'sQTL'
amd_joint_qqplot_data = rbind(amd_qqplot_data[str_detect(tissue,'(RPE|GWAS)'),],amd_sqtl_qqplot_data)
p3 = plot_joint_qqplot(amd_joint_qqplot_data)


myopia_sqtl_pval_list = make_sqtl_pval_list(myopia_gwas)
myopia_sqtl_qqplot_data = make_qqplot_data(myopia_sqtl_pval_list)
myopia_qqplot_data$type = 'eQTL'
myopia_sqtl_qqplot_data$type = 'sQTL'
myopia_joint_qqplot_data = rbind(myopia_qqplot_data[str_detect(tissue,'(RPE|GWAS)'),],myopia_sqtl_qqplot_data)
p4 = plot_joint_qqplot(myopia_joint_qqplot_data)

# make plot:
p = plot_grid(p3+ggtitle('AMD'),p4+ggtitle('Myopia'),nrow=1,labels=c('A','B'))
fig_fn = sprintf('%s/RPE_eQTL_sQTL_GWAS_qqplot.pdf',fig_dir)
save_plot(fig_fn,p,base_height=4,base_width=9)

fig_fn = sprintf('%s/RPE_eQTL_sQTL_GWAS_qqplot.png',fig_dir)
save_plot(fig_fn,p,base_height=4,base_width=9)


amd_ks_pval = foreach(pval = c(amd_pval_list[2:3],amd_sqtl_pval_list),.combine='c')%dopar%{
	result = ks.test(pval,amd_gwas$pval)
	result$p.value
}

amd_ks_pval # 0.0004956043 0.0090708669 0.0143265326 0.3113218811

myopia_ks_pval = foreach(pval = c(myopia_pval_list[2:3],myopia_sqtl_pval_list),.combine='c')%dopar%{
	result = ks.test(pval,myopia_gwas$pval)
	result$p.value
}

myopia_ks_pval # 0.02311454 0.29582928 0.53944768 0.58315872
