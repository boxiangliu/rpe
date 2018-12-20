library(data.table)
library(stringr)
library(cowplot)
library(ggplot2)
library(gap)
source('utils/genome_annotation.R')

fig_dir = '../figures/rasqual/quality_control/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

# Functions:
read_rasqual_chr22 = function(){
	tmp_file=tempfile()
	system(sprintf("cat ../processed_data/rasqual/output/glucose/joint/chr22/*.txt | awk '{print $1,$2,$3,$4,$10,$11,$12}' > %s",tmp_file))
	ras=fread(tmp_file,col.names=c('fid','sid','chr','pos','log10qval','chisq','pi'))
	file.remove(tmp_file)
	ras[,fid:=str_split_fixed(fid,'_',2)[,1]]
	ras[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
	return(ras)
}

plot_pval_histogram = function(ras){
	p = ggplot(ras,aes(pval)) + 
		geom_histogram(bins=1000,color='black',fill='black') + 
		xlab('p-value') + 
		ylab('Frequency')
	return(p)
}

plot_qqplot = function(ras){
	idx = sample(nrow(ras),size=1e5)
	plot.new()
	x = qqunif(ras$pval[idx],plot.it=FALSE)
	x = as.data.frame(x)
	p = ggplot(x,aes(x,y)) + 
		geom_point() + 
		geom_abline(intercept=0,slope=1,color='red')+
		xlab('Expected p-values') + 
		ylab('Observed p-values')
	return(p)
}

plot_pvalue_vs_position = function(ras){
	ras2 = ras[abs(dist)<1e6]
	set.seed(42)
	idx = sample(nrow(ras2),2e5)
	ras2 = ras2[idx,]
	p = ggplot(ras2,aes(dist,-log10(pval),size=(10*abs(pi-0.5))^2))+
		geom_point(alpha=0.2)+
		scale_size_continuous(name=expression('|'*pi*'-0.5|'),breaks=(10*c(0,0.1,0.2,0.3,0.4,0.5))^2,labels=c(0,0.1,0.2,0.3,0.4,0.5))+
		theme(legend.position=c(0.05,0.99),legend.justification=c('left','top'))+
		scale_y_sqrt()+
		xlab('Distance to TSS')+
		ylab(expression(-log[10]*'(p-value)'))
	return(p)
}

plot_eqtl_vs_dist = function(ras){
	sig_eqtl=ras[pval<1e-4,]
	p=ggplot(sig_eqtl,aes(dist))+
		geom_density(color='black',fill='black')+
		xlab('Distance to TSS')+
		ylab('Frequency')
	return(p)
}

# Read RASQUAL result: 
ras = read_rasqual_chr22()

# Plot histogram of p-value:
p1 = plot_pval_histogram(ras)
save_plot('../figures/rasqual/pvalue_distribution.pdf',p1)

# Plot QQ-plot:
p2 = plot_qqplot(ras)

# Plot p-value vs genomic position:
gencode = read_gencode(gencode_fn)
ras = merge(ras,gencode[,list(gene_id,start,end)],by.x='fid',by.y='gene_id')
ras[,dist:=pos-start]
p3 = plot_pvalue_vs_position(ras)

# Plot a histogram of eQTL along genomic position:
p4 = plot_eqtl_vs_dist(ras)


# Make entire plot:
top = plot_grid(p1,p2,labels=c('A','B'),nrow=1)
bottom = plot_grid(p3,p4,labels=c('C','D'),nrow=1)
p = plot_grid(top,bottom,nrow=2,labels='')
fig_fn = sprintf('%s/quality_control_uc.png',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)

top = plot_grid(p1,p2,labels=c('a','b'),nrow=1)
bottom = plot_grid(p3,p4,labels=c('c','d'),nrow=1)
p = plot_grid(top,bottom,nrow=2,labels='')
fig_fn = sprintf('%s/quality_control_lc.png',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)