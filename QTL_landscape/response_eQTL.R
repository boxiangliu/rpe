library(data.table)
library(stringr)
library(cowplot)
library(foreach)

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
expression_fn = '../data/rnaseq/count/merged/rpe.gene_count'
genotype_dir = '../data/genotype/filt/'
response_eQTL_fn = '../processed_data/QTL_landscape/eQTL_manhattan/response_eQTLs.txt'
glucose_covariate_fn = '../processed_data/select_covariates/merge_covariates/glu_sex_geno_sva.tsv'
galactose_covariate_fn = '../processed_data/select_covariates/merge_covariates/gal_sex_geno_sva.tsv'

fig_dir = '../figures/QTL_landscape/response_eQTL/'

if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_response_eQTL = function(fn){
	response_eQTL = fread(fn)
	return(response_eQTL)
}

find_rasqual_fn = function(gene_name,dir){
	fn = list.files(dir,pattern=paste0('_',gene_name,'.txt'),recursive=TRUE,full.names=TRUE)
	stopifnot(length(fn)==1)
	return(fn)
}

read_rasqual = function(fn){
	rasqual = fread(fn,header=FALSE)[,c(1,2,3,4,11)]
	colnames = c('fid','sid','chr','pos','chisq')
	setnames(rasqual,colnames)
	rasqual[,rank := rank(-chisq,ties.method='random')]
	rasqual = rasqual[rank == 1]
	rasqual$rank = NULL
	split_fid = str_split_fixed(rasqual$fid,'_',2)
	rasqual$gene_id = split_fid[,1]
	rasqual$gene_name = split_fid[,2]
	rasqual$fid = NULL
	rasqual[,pval := pchisq(chisq,df=1,lower.tail=FALSE)]
	rasqual[,logp := -log10(pval)]
	return(rasqual)
}


read_genotype = function(snp,dir = genotype_dir){
	chr = str_replace(snp$chr,'chr','')
	pos = snp$pos
	tmp_fn = tempfile()
	tmp_fn = paste0(tmp_fn,'.vcf')
	vcf_fn = sprintf('%s/rpe.imputed.chr%s.all_filters.vcf.gz',dir,chr)
	region = sprintf('%s:%s',chr,pos)
	command = sprintf('bcftools view -r %s %s > %s',region,vcf_fn,tmp_fn)
	system(command)
	command = sprintf('bcftools query -H -f "%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%DS]\\n" %s',tmp_fn)
	result = system(command,intern=TRUE)
	stopifnot(length(result)==2)
	split_header = str_split(result[1],'\t')[[1]]
	split_header = str_extract(split_header,'(?<=])(.+?)(?=(:|$))')
	split_result = str_split(result[2],'\t')[[1]]
	ref = split_result[3]
	alt = tolower(split_result[4])
	genotype = data.table(genotype = round(as.numeric(split_result[5:length(split_result)])))
	genotype[, genotype := ifelse(genotype == 0, paste0(ref,'/',ref), genotype)]
	genotype[, genotype := ifelse(genotype == 1, paste0(ref,'/',alt), genotype)]
	genotype[, genotype := ifelse(genotype == 2, paste0(alt,'/',alt), genotype)]
	genotype$sample = split_header[5:length(split_header)]
	unlink(tmp_fn)
	return(genotype)
}

update_name = function(x){
	dna2rna = fread('../data/meta/dna2rna.txt',colClasses = c('character','character'))
	idx = match(x$sample,dna2rna$DNA)
	x$sample = dna2rna[idx,RNA]
	return(x)
}


read_expression = function(fn){
	expression = fread(fn)
	return(expression)
}

subset_expression_matrix = function(expression_matrix,gene,treatment){
	expression = expression_matrix[str_detect(gene_id, gene),]
	expression = expression[,str_detect(colnames(expression),treatment),with=FALSE]
	colnames(expression) = str_replace(colnames(expression),treatment,'')
	colnames(expression) = str_replace(colnames(expression),'\\.$','')
	expression = data.table(sample=colnames(expression),expression=unlist(expression[1,]))
	return(expression)
}

read_covariate = function(fn){
	covariate = fread(fn)
	return(covariate)
}

get_residual = function(expression,covariate){
	design = merge(expression,covariate)
	sample = design$sample
	design$sample = NULL
	result = lm(expression ~ .,data = design)
	residual = result$residuals
	residual = data.table(sample, residual)
	return(residual)
}

read_ase = function(){
	return(ase)
}

plot_eQTL = function(gxe){
	p = ggplot(gxe,aes(genotype,residual,fill = treatment)) + 
		geom_boxplot() +
		facet_grid(.~treatment) +
		xlab('Genotype') + 
		ylab('Expression residuals') + 
		scale_fill_discrete(guide='none') + 
		theme(strip.background =element_rect(fill="white",color='black',size=1,linetype = 'solid'))
	return(p)
}

#--------- Main ------------# 
response_eQTL = read_response_eQTL(response_eQTL_fn)
glucose_eQTL = response_eQTL[eQTL == 'Glucose-specific']
galactose_eQTL = response_eQTL[eQTL == 'Galactose-specific']

glucose_covariate = read_covariate(glucose_covariate_fn)
galactose_covariate = read_covariate(galactose_covariate_fn)

galactose_p = foreach(i = 1:10)%do%{

	gene_id = galactose_eQTL[i,gene_id]
	gene_name = galactose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	fn = find_rasqual_fn(gene_name,glucose_dir)
	rasqual = read_rasqual(fn)
	snp = rasqual[,list(sid,chr,pos)]

	genotype = read_genotype(snp)
	genotype = update_name(genotype)

	expression_matrix = read_expression(expression_fn)
	glucose_expression = subset_expression_matrix(expression_matrix,gene_id,'glucose')
	glucose_residual = get_residual(glucose_expression,glucose_covariate)
	
	galactose_expression = subset_expression_matrix(expression_matrix,gene_id,'galactose')
	galactose_residual = get_residual(galactose_expression,galactose_covariate)

	glucose_residual$treatment = 'Glucose'
	galactose_residual$treatment = 'Galactose'
	residual = rbind(glucose_residual,galactose_residual)
	gxe = merge(genotype,residual)
	p = plot_eQTL(gxe)

	fig_fn = sprintf('%s/galactoseSpecific-%s-%s-%s.pdf',fig_dir,gene_name,gene_id,snp$sid) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}

glucose_p = foreach(i = 1:10)%do%{

	gene_id = glucose_eQTL[i,gene_id]
	gene_name = glucose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	fn = find_rasqual_fn(gene_name,glucose_dir)
	rasqual = read_rasqual(fn)
	snp = rasqual[,list(sid,chr,pos)]

	genotype = read_genotype(snp)
	genotype = update_name(genotype)

	expression_matrix = read_expression(expression_fn)
	glucose_expression = subset_expression_matrix(expression_matrix,gene_id,'glucose')
	glucose_residual = get_residual(glucose_expression,glucose_covariate)
	
	galactose_expression = subset_expression_matrix(expression_matrix,gene_id,'galactose')
	galactose_residual = get_residual(galactose_expression,galactose_covariate)

	glucose_residual$treatment = 'Glucose'
	galactose_residual$treatment = 'Galactose'
	residual = rbind(glucose_residual,galactose_residual)
	gxe = merge(genotype,residual)
	p = plot_eQTL(gxe)

	fig_fn = sprintf('%s/glucoseSpecific-%s-%s-%s.pdf',fig_dir,gene_name,gene_id,snp$sid) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}