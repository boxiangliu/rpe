library(data.table)
library(stringr)
library(foreach)
library(cowplot)
source('utils/gtex_tissue_info.R')

# Make QQ plot:
gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'
specific_gene_dir = '../processed_data/rpe_specific_genes.GTExV7/'
fn_list = list.files(specific_gene_dir,'all_genes',full.names=TRUE)
fig_dir = '../figures/rpe_specific_genes/plot_QQplot/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

read_zscore = function(fn_list){
	zscore = foreach(fn = fn_list,.combine='rbind') %dopar%{
		tissue = str_split_fixed(basename(fn),'\\.',3)[,2]
		tissue = str_replace_all(tissue,'_',' ')
		tissue_id = tissue_to_tissue_id[tissue]
		x = fread(fn)
		x$tissue = tissue_id
		return(x)
	}
	return(zscore)
}

zscore = read_zscore(fn_list)

get_mean_zscore = function(zscore){
	set.seed(42)
	zscore[,rank:=rank(zscore,ties.method='random'),by='tissue']
	mean_zscore = zscore[,list(mean_zscore = mean(zscore,na.rm=TRUE)),by='rank']
	return(mean_zscore)
}

mean_zscore = get_mean_zscore(zscore)

make_qqplot_data = function(zscore,mean_zscore){
	qqplot_data = foreach(t = unique(zscore$tissue),.combine='rbind')%do%{
		x = qqplot(mean_zscore$mean_zscore,zscore[tissue==t,zscore],plot.it=FALSE)
		x = data.frame(x)
		x$tissue = t
		return(x)
	}
	setDT(qqplot_data)
	qqplot_data[,alpha := ifelse(tissue=='RPE','high','low')]
	return(qqplot_data)
}

qqplot_data = make_qqplot_data(zscore,mean_zscore)

plot_qqplot = function(qqplot_data){
	p = ggplot(qqplot_data_small,aes(x,y,color=tissue,alpha=alpha))+
		geom_abline(intercept=0,slope=1,color='black',linetype='dashed') + 
		geom_line() +
		scale_color_manual(values=tissue_color,name='',labels=tissue_abbreviation) + 
		scale_alpha_manual(values=c(high=1,low=0.2),guide='none') + 
		theme(legend.position='bottom',legend.text=element_text(size=6)) + 
		xlab('Expected z-score') + 
		ylab('Observed z-score') + 
		annotate(geom = 'text',x=0,y=3.5,label='Testis')
	return(p)
}

p = plot_qqplot(qqplot_data)
fig_fn = sprintf('%s/qqplot.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=6)
