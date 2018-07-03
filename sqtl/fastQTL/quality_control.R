library(data.table)
library(stringr)
library(cowplot)
library(gap)

# Variables:
glucose_fn = '../processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz'
fig_dir = '../figures/sqtl/fastQTL/quality_control/'
if (!dir.exists(fig_dir)){dir.create(fig_dir,recursive=TRUE)}

# Functions:
read_sqtl = function(fn){
	sqtl=fread(sprintf('zcat %s',fn),header=F)
	setnames(sqtl,c('intron','snp','dist','pval','beta','varbeta'))
	return(sqtl)
}


plot_pval_histogram = function(fas){
	idx = sample(nrow(fas),size=1e5)
	p = ggplot(fas[idx,],aes(pval)) + 
		geom_histogram(bins=1000,color='black',fill='black')+
		xlab('p-value')+
		ylab('Frequency')
	return(p)
}

plot_qqplot = function(fas){
	idx = sample(nrow(fas),size=1e5)
	plot.new()
	x = qqunif(fas$pval[idx],plot.it=FALSE)
	x = as.data.frame(x)
	p = ggplot(x,aes(x,y)) + 
		geom_point() + 
		geom_abline(intercept=0,slope=1,color='red')+
		xlab('Expected p-values') + 
		ylab('Observed p-values')
	return(p)
}


plot_sqtl_vs_dist = function(sqtl){
	sig_sqtl=sqtl[pval<1e-4,]
	p=ggplot(sig_sqtl,aes(dist))+
		geom_histogram(bins=100,color='black',fill='black')+
		xlab('Distance to splicing donor/acceptor')+
		ylab('Frequency')
	return(p)
}


plot_sqtl_vs_intron_boundary = function(sqtl){
	sig_sqtl=sqtl[pval<1e-4,]
	temp=str_split_fixed(sig_sqtl$intron,":",4)
	sig_sqtl$start=as.numeric(temp[,2])
	sig_sqtl$end=as.numeric(temp[,3])
	sig_sqtl=sig_sqtl[,length:=end-start]
	sig_sqtl=sig_sqtl[,fraction:=dist/length]
	sig_sqtl_within_intron=sig_sqtl[fraction<=1&fraction>=0,]
	p=ggplot(sig_sqtl_within_intron,aes(fraction))+
	geom_histogram(bins=100,color='black',fill='black')+
	xlab('Fraction distance within intron') + 
	ylab('Frequency')
	return(p)
}

# Main:
glucose = read_sqtl(glucose_fn)
p1 = plot_pval_histogram(glucose)
p2 = plot_qqplot(glucose)
p3 = plot_sqtl_vs_dist(glucose)
p4 = plot_sqtl_vs_intron_boundary(glucose)

p = plot_grid(p1,p2,p3,p4,labels=LETTERS[1:4],nrow=2)
fig_fn = sprintf('%s/glucose_sqtl_quality_control_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)
fig_fn = sprintf('%s/glucose_sqtl_quality_control_uc.png',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)

p = plot_grid(p1,p2,p3,p4,labels=letters[1:4],nrow=2)
fig_fn = sprintf('%s/glucose_sqtl_quality_control_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)
fig_fn = sprintf('%s/glucose_sqtl_quality_control_lc.png',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)