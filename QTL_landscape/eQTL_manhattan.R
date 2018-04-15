library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)
library(manhattan)
library(cowplot)

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
treeQTL_MT_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
fig_dir = '../figures/QTL_landscape/eQTL_manhattan/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_rasqual = function(fn){
	rasqual = fread(fn,header=FALSE)[,c(1,3,4,11)]
	colnames = c('fid','chr','pos','chisq')
	setnames(rasqual,colnames)
	rasqual[,rank := rank(-chisq,ties.method='random')]
	rasqual = rasqual[rank == 1]
	rasqual$rank = NULL
	rasqual[,gene_id := str_split_fixed(fid,'_',2)[,1]]
	rasqual$fid = NULL
	rasqual[,pval := pchisq(chisq,df=1,lower.tail=FALSE)]
	rasqual[,logp := -log10(pval)]
	return(rasqual)
}

read_TreeQTL_MT = function(fn){
	treeQTL_MT = fread(fn)
	return(treeQTL_MT)
}

label_eQTL = function(gene_id, eqtl_list){
	eqtl_indicator = gene_id %in% eqtl_list
	return(eqtl_indicator)
}

rbind_eQTLs = function(x,y){
	z = copy(y)
	z[,logp := -logp]
	rbind(x,z)
}

plot_manhattan = function(eqtl){
	data = copy(eqtl)
	build = 'hg19'
	data[,y:=logp]
	data[,chrom:=chr]
	data = add_cumulative_pos(data, build)
	color1 = 'black'
	color2 = 'grey'
	data = add_color(data, color1 = color1, color2 = color2)
	data[eQTL==TRUE,color := 'red']
	data[treatment_specific == TRUE,color:='blue']
	chrom_lengths = get_chrom_lengths(build)
	xmax = get_total_length(chrom_lengths)
	x_breaks = get_x_breaks(chrom_lengths)
	color_map = unique(data$color)
	names(color_map) = unique(data$color)

	p=ggplot(data, aes(x = cumulative_pos, y = y, color = color)) + 
		geom_point(alpha = 0.5) + 
		theme_classic() + 
		scale_x_continuous(limits = c(0, xmax), expand = c(0.01, 0), breaks = x_breaks, labels = names(x_breaks), name = "Chromosome") + 
		scale_y_continuous(expand = c(0.01, 0), name = expression("-log10(P-value)"),labels = abs) + 
		scale_color_manual(values = color_map, name = '', breaks = c('red','blue'), labels = c('eQTL','Treatment-specific eQTL')) + 
		theme(legend.position = 'top') +
		facet_grid(treatment~., scales = 'free_y') 
	return(p)
}

plot_scatterplot = function(glucose,galactose){
	data = merge(
		glucose[,list(gene_id,glucose_logp=logp,glucose_specific=treatment_specific)],
		galactose[,list(gene_id,galactose_logp=logp,galactose_specific=treatment_specific)]
		,by='gene_id')
	data$color = 'black'
	data[,color := ifelse(glucose_specific,'red',color)]
	data[,color := ifelse(galactose_specific,'blue',color)]
	p = ggplot(data,aes(glucose_logp,galactose_logp,color=color)) + 
		geom_point(alpha = 0.5) + 
		geom_abline(intercept = 0, slope =1, color = 'black', linetype = 'dashed') + 
		scale_color_manual(name = '', breaks = c('blue','red'), labels = c('Galactose-specific','Glucose-specific'),values = c(black='black',red='red',blue='blue')) + 
		scale_x_log10() + 
		scale_y_log10() + 
		xlab(expression(Glucose -log['10'](P-value))) + 
		ylab(expression(Galactose -log['10'](P-value))) + 
		theme(legend.position = 'top') + 

	return(p)
}

# Main:
treeQTL_MT = read_TreeQTL_MT(treeQTL_MT_fn)

glucose_list = list.files(glucose_dir,pattern = 'txt',recursive=TRUE,full.names=TRUE)
glucose = foreach(fn = glucose_list, .combine = rbind)%dopar%{
	rasqual = read_rasqual(fn)
	return(rasqual)
}

glucose$eQTL = label_eQTL(glucose$gene_id,treeQTL_MT[glucose==1,gene])
glucose$treatment_specific = label_eQTL(glucose$gene_id,treeQTL_MT[glucose==1&galactose==0,gene])
glucose$treatment = 'Glucose'

galactose_list = list.files(galactose_dir,pattern = 'txt',recursive=TRUE,full.names=TRUE)
galactose = foreach(fn = galactose_list, .combine = rbind)%dopar%{
	rasqual = read_rasqual(fn)
	return(rasqual)
}
galactose$eQTL = label_eQTL(galactose$gene_id,treeQTL_MT[galactose==1,gene])
galactose$treatment_specific = label_eQTL(galactose$gene_id,treeQTL_MT[galactose==1&glucose==0,gene])
galactose$treatment = 'Galactose'

eqtl = rbind_eQTLs(galactose,glucose)
p = plot_manhattan(eqtl)
fig_fn = sprintf('%s/manhattan.pdf',fig_dir)
pdf(fig_fn, width = 8, height = 4)
print(p)
dev.off()

p = plot_scatterplot(glucose,galactose)
fig_fn = sprintf('%s/scatterplot.pdf',fig_dir)
pdf(fig_fn,width = 6, height = 6)
print(p)
dev.off()