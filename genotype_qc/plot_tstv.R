library(data.table)
library(foreach)
library(cowplot)

fig_dir = '../figures/genotype_qc/plot_tstv/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_tstv = function(fn){
	tstv = read.table(fn,nrows=1,skip=19,header=FALSE)
	setDT(tstv)
	tstv = tstv[,3:5]
	setnames(tstv,c('ts','tv','ts/tv'))
	return(tstv)
}


orig = foreach(i = 1:22, .combine = rbind)%do%{
	fn = sprintf('../data/genotype/orig/qc/rpe.imputed.chr%s.vcf.gz.stats',i)
	tstv = read_tstv(fn)
	tstv$chr = i
	return(tstv)
}
orig$type = 'orig'

filt = foreach(i = 1:22, .combine = rbind)%do%{
	fn = sprintf('../data/genotype/filt/qc/rpe.imputed.chr%s.all_filters.vcf.gz.stats',i)
	tstv = read_tstv(fn)
	tstv$chr = i
	return(tstv)
}
filt$type = 'filt'
tstv = rbind(orig,filt)
tstv[,chr := factor(chr,1:22)]

plot_tstv = function(tstv){
	p = ggplot(tstv,aes(chr,`ts/tv`,fill=type))+
		geom_bar(stat='identity',position='dodge') + 
		geom_hline(yintercept = 2.0, color = 'red', linetype = 'dashed') + 
		geom_hline(yintercept = 2.1, color = 'red', linetype = 'dashed') + 
		scale_y_continuous(breaks = c(seq(0,2.5,0.5),2.1), labels = c(seq(0,2.5,0.5),2.1)) + 
		scale_fill_discrete(name = 'VCF', breaks = c('filt','orig'), labels = c('Filtered','Original'))+
		xlab('Chromosome') + 
		ylab('Ts/Tv ratio')
	return(p)
}

p = plot_tstv(tstv)
fig_fn = sprintf('%s/tstv.pdf',fig_dir)
save_plot(fig_fn, p, base_width = 7, base_height = 5)

tstv[chr == 8 & type == 'filt'] # 1.98
tstv[chr == 16 & type == 'filt'] # 1.93
