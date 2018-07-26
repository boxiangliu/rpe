library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)
source('utils/plot_ase.R')

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
expression_fn = '../data/rnaseq/count/merged/rpe.gene_count'
genotype_dir = '../data/genotype/filt/'
glucose_eQTL_fn = '../processed_data/QTL_landscape/eQTL_manhattan/rasqual_glucose_top_eQTL.txt'
galactose_eQTL_fn = '../processed_data/QTL_landscape/eQTL_manhattan/rasqual_galactose_top_eQTL.txt'
glucose_covariate_fn = '../processed_data/select_covariates/merge_covariates/glu_sex_geno_sva.tsv'
galactose_covariate_fn = '../processed_data/select_covariates/merge_covariates/gal_sex_geno_sva.tsv'
imprinted_gene_fn = '../data/imprinted_genes/imprinted_genes.txt'
fig_dir = '../figures/QTL_landscape/response_eQTL/'
out_dir = '../processed_data/QTL_landscape/response_eQTL/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


read_response_eQTL = function(fn){
	response_eQTL = fread(fn)
	return(response_eQTL)
}

read_genotype = function(snp,dir = genotype_dir){
	chr = str_replace(snp$chr,'chr','')
	pos = snp$pos
	ref = snp$ref
	alt = snp$alt

	tmp_fn = tempfile()
	tmp_fn = paste0(tmp_fn,'.vcf')
	vcf_fn = sprintf('%s/rpe.imputed.chr%s.all_filters.vcf.gz',dir,chr)
	region = sprintf('%s:%s',chr,pos)
	command = sprintf('bcftools view -r %s %s > %s',region,vcf_fn,tmp_fn)
	system(command)

	command = sprintf('bcftools query -H -f "%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%DS]\\n" %s',tmp_fn)
	result = system(command,intern=TRUE)
	split_header = str_split(result[1],'\t')[[1]]
	split_header = str_extract(split_header,'(?<=])(.+?)(?=(:|$))')
	pattern = paste0(paste(chr,pos,ref,alt,sep='\t'),'\t')
	idx = which(str_detect(result,pattern))
	stopifnot(length(idx)==1)
	split_result = str_split(result[idx],'\t')[[1]]

	alt = tolower(alt)
	genotype = data.table(genotype = round(as.numeric(split_result[5:length(split_result)])))
	g0 = paste0(ref,ref)
	g1 = paste0(ref,alt)
	g2 = paste0(alt,alt)
	genotype[, genotype := ifelse(genotype == 0, g0, genotype)]
	genotype[, genotype := ifelse(genotype == 1, g1, genotype)]
	genotype[, genotype := ifelse(genotype == 2, g2, genotype)]
	genotype[, genotype := factor(genotype,levels = c(g0,g1,g2))]
	genotype$sample = split_header[5:length(split_header)]
	unlink(tmp_fn)
	return(genotype)
}

update_name = function(x){
	dna2rna = fread('../data/meta/dna2rna.txt',colClasses = c('character','character','character'))
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

normalize_residual = function(residual){
	max = max(residual$residual)
	min = min(residual$residual)
	range = max-min
	residual[,residual:=residual/range]
	return(residual)
}

plot_eQTL = function(gxe){
	p = ggplot(gxe,aes(genotype,residual,fill = treatment)) + 
		geom_boxplot() +
		facet_grid(.~treatment) +
		xlab('') + 
		ylab('Normalized residuals') + 
		scale_fill_manual(guide='none',values = c(Glucose = 'red', Galactose = 'green')) + 
		theme(strip.background =element_rect(fill="white",color='black',size=1,linetype = 'solid'))
	return(p)
}


get_snp = function(rasqual){
	snp = rasqual[,list(sid,chr,pos)]
	split_sid = str_split_fixed(snp$sid,'_',5)
	snp$ref = split_sid[,3]
	snp$alt = split_sid[,4]
	return(snp)
}

plot_eQTL_driver = function(gene_id,gene_name,dir,expression_matrix){
	fn = find_rasqual_fn(gene_id,gene_name,dir)
	rasqual = read_rasqual(fn)
	snp = get_snp(rasqual)

	genotype = read_genotype(snp)
	genotype = update_name(genotype)

	glucose_expression = subset_expression_matrix(expression_matrix,gene_id,'glucose')
	glucose_residual = get_residual(glucose_expression,glucose_covariate)
	glucose_residual = normalize_residual(glucose_residual)

	galactose_expression = subset_expression_matrix(expression_matrix,gene_id,'galactose')
	galactose_residual = get_residual(galactose_expression,galactose_covariate)
	galactose_residual = normalize_residual(galactose_residual)

	glucose_residual$treatment = 'Glucose'
	galactose_residual$treatment = 'Galactose'
	residual = rbind(glucose_residual,galactose_residual)
	gxe = merge(genotype,residual)
	gxe[,treatment:=factor(treatment,c('Glucose','Galactose'))]
	p = plot_eQTL(gxe)
	return(p)
}

merge_glucose_galactose = function(glucose,galactose){
	data = merge(
		glucose[,list(gene_id,gene_name,glucose_pi = pi, glucose_logp=logp,glucose_specific=treatment_specific,glucose_shared = eQTL&!treatment_specific)],
		galactose[,list(gene_id,galactose_pi = pi, galactose_logp=logp,galactose_specific=treatment_specific,galactose_shared = eQTL&!treatment_specific)]
		,by='gene_id')
	return(data)
}

rank_eQTL_by_pi = function(data){
	delta_pi = abs(data$glucose_pi - 0.5) - abs(data$galactose_pi - 0.5)
	delta_pi = ifelse(data$glucose_specific==TRUE,delta_pi,-Inf)
	data$glucose_specific_rank = rank(-delta_pi,ties.method='first')

	delta_pi = abs(data$galactose_pi - 0.5) - abs(data$glucose_pi - 0.5)
	delta_pi = ifelse(data$galactose_specific==TRUE,delta_pi,-Inf)
	data$galactose_specific_rank = rank(-delta_pi,ties.method='first')

	sum_pi = abs(data$glucose_pi - 0.5) + abs(data$galactose_pi - 0.5)
	sum_pi = ifelse(data$glucose_shared==TRUE,sum_pi, -Inf)
	data$shared_rank = rank(-sum_pi,ties.method='first')
	return(data)
}

#--------- Main ------------# 
glucose = fread(glucose_eQTL_fn)
galactose = fread(galactose_eQTL_fn)
glucose_eQTL = glucose[treatment_specific==TRUE]
galactose_eQTL = galactose[treatment_specific==TRUE]
shared_eQTL = glucose[eQTL==TRUE&treatment_specific==FALSE]
glucose_covariate = read_covariate(glucose_covariate_fn)
galactose_covariate = read_covariate(galactose_covariate_fn)
expression_matrix = read_expression(expression_fn)

#--------------------#
# galactose-specific #
#--------------------#
galactose_p = foreach(i = 1:nrow(galactose_eQTL))%do%{

	gene_id = galactose_eQTL[i,gene_id]
	gene_name = galactose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,galactose_dir,expression_matrix) + ggtitle(gene_name)
	fig_fn = sprintf('%s/galactose_specific/galactoseSpecific-%s-%s.pdf',fig_dir,gene_name,gene_id) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}
names(galactose_p) = galactose_eQTL$gene_name
out_fn = sprintf('%s/galactose_specific_eQTL_plots.rds',out_dir)
saveRDS(galactose_p,out_fn)

#------------------#
# glucose-specific #
#------------------#
glucose_p = foreach(i = 1:nrow(glucose_eQTL))%do%{

	gene_id = glucose_eQTL[i,gene_id]
	gene_name = glucose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,glucose_dir,expression_matrix) + ggtitle(gene_name)
	fig_fn = sprintf('%s/glucose_specific/glucoseSpecific-%s-%s.pdf',fig_dir,gene_name,gene_id) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}
names(glucose_p) = glucose_eQTL$gene_name
out_fn = sprintf('%s/glucose_specific_eQTL_plots.rds',out_dir)
saveRDS(glucose_p,out_fn)

#--------#
# shared #
#--------#
shared_p = foreach(i = 1:nrow(shared_eQTL))%do%{

	gene_id = shared_eQTL[i,gene_id]
	gene_name = shared_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,glucose_dir,expression_matrix) + ggtitle(gene_name)
	fig_fn = sprintf('%s/shared/shared-%s-%s.pdf',fig_dir,gene_name,gene_id) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}
names(shared_p) = shared_eQTL$gene_name
out_fn = sprintf('%s/shared_eQTL_plots.rds',out_dir)
saveRDS(shared_p,out_fn)

#----------#
# ASE plot #
#----------#
imprinted_gene = fread(imprinted_gene_fn)
data = merge_glucose_galactose(glucose,galactose)
data = data[!(gene_name %in% imprinted_gene$Gene)]
data = rank_eQTL_by_pi(data)
top_tissue_specific_eqtl = data[glucose_specific_rank<=5]
setorder(top_tissue_specific_eqtl,glucose_specific_rank)
plot_data=get_ase_plot_data(top_tissue_specific_eqtl)
glucose_ase_plot = plot_ase(plot_data) 

top_tissue_specific_eqtl = data[galactose_specific_rank<=5]
setorder(top_tissue_specific_eqtl,glucose_specific_rank)
plot_data=get_ase_plot_data(top_tissue_specific_eqtl)
galactose_ase_plot = plot_ase(plot_data)

top_shared_eqtl = data[shared_rank<=5]
setorder(top_shared_eqtl,shared_rank)
plot_data=get_ase_plot_data(top_shared_eqtl)
shared_ase_plot = plot_ase(plot_data)

p = plot_grid(galactose_ase_plot,glucose_ase_plot,shared_ase_plot,nrow=3,align='v',labels=c('A','B','C'))
fig_fn = sprintf('%s/ase/ase_plot_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_height = 8, base_width = 8)

p = plot_grid(galactose_ase_plot,glucose_ase_plot,shared_ase_plot,nrow=3,align='v',labels=c('a','b','c'))
fig_fn = sprintf('%s/ase/ase_plot_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_height = 8, base_width = 8)