library(data.table)
library(foreach)
library(cowplot)

glucose_dir='../processed_data/sqtl/optimal_covariate_sva/test_covariate/glucose/'
galactose_dir='../processed_data/sqtl/optimal_covariate_sva/test_covariate/galactose/'
fig_dir = '../figures/sqtl/optimal_covariate_sva/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

count_significant_sqtl = function(in_dir){
	sig = foreach(i = 0:6,.combine='rbind')%do%{
		fn = sprintf('%s/chr1.permute.%s.txt.gz',in_dir,i)
		sqtl=fread(sprintf('zcat %s',fn),select=c(1,6,7,16),col.names=c('cluster','snp','dist','pval'))
		sqtl[,fdr := p.adjust(pval)]
		fdr5 = sqtl[,sum(fdr<0.05,na.rm=TRUE)]
		fdr10 = sqtl[,sum(fdr<0.10,na.rm=TRUE)]
		x = data.table(sig = c(fdr5,fdr10), fdr = c(0.05,0.1), cov = i)
		return(x)
	}
	sig[,fdr:=as.factor(fdr)]
	return(sig)
}

plot_siginificant_sqtl = function(sig,title){
	p = ggplot(sig,aes(cov,sig,color=fdr)) + 
		geom_point() + 
		geom_line() + 
		scale_x_continuous(breaks=0:6,labels=c('intercept','sex',paste0('geno PC',1:3),paste0('SV',1:2)))+
		theme(axis.text.x=element_text(angle=45,hjust=1),legend.position=c(0.95,0.95),legend.justification=c('right','top'))+
		scale_color_discrete(name='FDR')+
		xlab('') + 
		ylab('# significant sQTL clusters\n(FDR < 0.05 or 0.1)') + 
		ggtitle(title)
	return(p)
}

sig = count_significant_sqtl(glucose_dir)
p1 = plot_siginificant_sqtl(sig,'Glucose')

sig = count_significant_sqtl(galactose_dir)
p2 = plot_siginificant_sqtl(sig,'Galactose')

p = plot_grid(p1,p2,labels=c('A','B'))
save_plot(sprintf('%s/sqtl_vs_cov.pdf',fig_dir),p,base_width=8)