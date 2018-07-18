library(data.table)
library(stringr)
library(cowplot)

intron_count_fn = '../processed_data/sqtl/visualization/prepare_results/rdh5_perind_numers.counts.gz'
genotype_dir = '../data/genotype/filt/'
expression_fn = '../data/rnaseq/count/merged/rpe.gene_count'
glucose_covariate_fn = '../processed_data/select_covariates/merge_covariates/glu_sex_geno_sva.tsv'
fig_dir = '../figures/figure5/qtl/'
out_dir = '../processed_data/figure5/qtl/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

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
	mean = mean(residual$residual)
	range = max-min
	residual[,residual:=(residual-mean)/range]
	return(residual)
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

get_splicing_ratio = function(){
	intron_count = read.table(intron_count_fn,check.names=FALSE)
	splicing_ratio = intron_count['chr12:56115278:56117670:clu_1202',]/colSums(intron_count)
	splicing_ratio = data.frame(t(splicing_ratio),check.names=FALSE)
	splicing_ratio$sample = rownames(splicing_ratio)
	setnames(splicing_ratio,'chr12:56115278:56117670:clu_1202','residual')
	setDT(splicing_ratio)
	return(splicing_ratio)
}

get_genotype = function(){
	snp=list(chr='chr12',pos=56115778,ref='C',alt='A')
	genotype = read_genotype(snp)
	genotype = update_name(genotype)
	return(genotype)
}

get_expression_residual = function(){
	expression_matrix = read_expression(expression_fn)
	glucose_expression = subset_expression_matrix(expression_matrix,'ENSG00000135437','glucose')
	glucose_covariate = read_covariate(glucose_covariate_fn)
	glucose_residual = get_residual(glucose_expression,glucose_covariate)
	glucose_residual = normalize_residual(glucose_residual)
	return(glucose_residual)
}


get_plot_data = function(){
	to_plot1 = merge(genotype,splicing_ratio,by='sample')
	to_plot1$type='sQTL'
	to_plot2 = merge(genotype,glucose_residual,by='sample')
	to_plot2$type='eQTL'
	to_plot = rbind(to_plot1,to_plot2)
	to_plot[,type:=factor(type,levels=c('sQTL','eQTL'))]
	return(to_plot)
}

plot_boxplot = function(){
	p = ggplot(to_plot,aes(genotype,residual)) + 
		geom_boxplot() + 
		facet_wrap(~type,scales='free') + 
		xlab(NULL) + 
		ylab(NULL) + 
		theme(strip.background = element_rect(colour = "black", fill = "white"))
	return(p)
}

splicing_ratio = get_splicing_ratio()
genotype = get_genotype()
glucose_residual = get_expression_residual()
to_plot = get_plot_data()
p = plot_boxplot()

fig_fn = sprintf('%s/qtl.pdf',fig_dir)
save_plot(fig_fn,p,base_height=2,base_width=2.2)
out_fn = sprintf('%s/qtl.rda',out_dir)
saveRDS(p,out_fn)