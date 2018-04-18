library(data.table)
library(stringr)
library(cowplot)
library(foreach)
library(doMC)
registerDoMC(15)

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
	rasqual = fread(fn,header=FALSE)[,c(1,2,3,4,11,12)]
	colnames = c('fid','sid','chr','pos','chisq','pi')
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
	g0 = paste0(ref,'/',ref)
	g1 = paste0(ref,'/',alt)
	g2 = paste0(alt,'/',alt)
	genotype[, genotype := ifelse(genotype == 0, g0, genotype)]
	genotype[, genotype := ifelse(genotype == 1, g1, genotype)]
	genotype[, genotype := ifelse(genotype == 2, g2, genotype)]
	genotype[, genotype := factor(genotype,levels = c(g0,g1,g2))]
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



plot_eQTL_driver = function(gene_id,gene_name,dir,expression_matrix){
	fn = find_rasqual_fn(gene_name,dir)
	rasqual = read_rasqual(fn)
	snp = rasqual[,list(sid,chr,pos)]

	genotype = read_genotype(snp)
	genotype = update_name(genotype)

	glucose_expression = subset_expression_matrix(expression_matrix,gene_id,'glucose')
	glucose_residual = get_residual(glucose_expression,glucose_covariate)
	
	galactose_expression = subset_expression_matrix(expression_matrix,gene_id,'galactose')
	galactose_residual = get_residual(galactose_expression,galactose_covariate)

	glucose_residual$treatment = 'Glucose'
	galactose_residual$treatment = 'Galactose'
	residual = rbind(glucose_residual,galactose_residual)
	gxe = merge(genotype,residual)
	p = plot_eQTL(gxe)
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

extract_ase = function(rasqual){
	ase = rasqual$pi
	pval = rasqual$pval
	zscore = abs(qnorm(pval/2))
	se = abs(ase-0.5)/zscore
	ci = se*1.96
	return(list(ase = ase,se = se,ci = ci))
}

read_ase = function(){
	return(ase)
}

make_ase_dataframe = function(glucose_ase,galactose_ase){
	glucose = data.table(treatment = 'Glucose',ase = glucose_ase$ase, ci = glucose_ase$ci)
	galactose = data.table(treatment = 'Galactose', ase = galactose_ase$ase, ci = galactose_ase$ci)
	ase = rbind(glucose,galactose)
	return(ase)
}

get_ase_plot_data = function(top_tissue_specific_eqtl){
	plot_data=foreach(i = 1:nrow(top_tissue_specific_eqtl),.combine = 'rbind') %do%{
		gene_id = top_tissue_specific_eqtl[i,gene_id]
		gene_name = top_tissue_specific_eqtl[i,gene_name]

		message(gene_id)
		message(gene_name)

		fn = find_rasqual_fn(gene_name,glucose_dir)
		rasqual = read_rasqual(fn)
		glucose_ase = extract_ase(rasqual)

		fn = find_rasqual_fn(gene_name,galactose_dir)
		rasqual = read_rasqual(fn)
		galactose_ase = extract_ase(rasqual)

		ase = make_ase_dataframe(glucose_ase,galactose_ase)
		ase$gene_name = gene_name
		return(ase)
	}
	plot_data[,gene_name := factor(gene_name,unique(gene_name))]
	return(plot_data)
}

plot_ase = function(ase){
	ggplot(ase,aes(treatment,ase,ymin = ase-ci,ymax = ase+ci, color = treatment)) + 
		geom_pointrange() + 
		geom_hline(yintercept = 0.5, color = 'red', linetype = 'dashed') + 
		xlab('') + 
		ylab('Allelic ratio (95% CI)') + 
		scale_color_discrete(guide = 'none') +
		facet_grid(.~gene_name) + 
		theme_bw() + 
		theme(strip.background=element_rect(color = 'black', fill = 'white', linetype = 'solid', size = 1))
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


galactose_p = foreach(i = 1:nrow(galactose_eQTL))%dopar%{

	gene_id = galactose_eQTL[i,gene_id]
	gene_name = galactose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,glucose_dir,expression_matrix)
	fig_fn = sprintf('%s/galactose_specific/galactoseSpecific-%s-%s-%s.pdf',fig_dir,gene_name,gene_id,snp$sid) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}

#------------------#
# glucose-specific #
#------------------#
glucose_p = foreach(i = 1:nrow(glucose_eQTL))%dopar%{

	gene_id = glucose_eQTL[i,gene_id]
	gene_name = glucose_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,galactose_dir,expression_matrix)
	fig_fn = sprintf('%s/glucose_specific/glucoseSpecific-%s-%s-%s.pdf',fig_dir,gene_name,gene_id,snp$sid) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}

#--------#
# shared #
#--------#
shared_p = foreach(i = 1:nrow(shared_eQTL))%dopar%{

	gene_id = shared_eQTL[i,gene_id]
	gene_name = shared_eQTL[i,gene_name]

	message(gene_id)
	message(gene_name)

	p = plot_eQTL_driver(gene_id,gene_name,glucose_dir,expression_matrix)
	fig_fn = sprintf('%s/shared/shared-%s-%s-%s.pdf',fig_dir,gene_name,gene_id,snp$sid) 
	save_plot(fig_fn,p,base_height=4,base_width=4)

	return(p)
}

#######################
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
fig_fn = sprintf('%s/ase/ase_plot.pdf',fig_dir)
save_plot(fig_fn,p,base_height = 8, base_width = 8)