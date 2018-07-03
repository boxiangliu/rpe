library(data.table)
library(foreach)
library(stringr)
library(cowplot)

star_alignment_dir = '../data/rnaseq/star_qc/'
fig_dir = '../figures/rnaseq_qc/plot_mapping_stats/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_STAR_mapping_stats = function(fn){
	lines = readLines(fn)
	stats = foreach(line = lines,.combine=rbind)%do%{
		if (str_detect(line,'\\|')){
			split_line = str_split_fixed(line,'\\|', 2)
			field = trimws(split_line[,1])
			value = trimws(split_line[,2])
			return(data.table(field,value))
		} else {
			return(NULL)
		}
	}
	return(stats)
}

plot_STAR_mapping_stats = function(x,xlabel){
	x = copy(x)
	setorder(x,value)
	x[,sample := factor(sample,unique(sample))]
	p = ggplot(x, aes(x = sample, y = value, fill = condition)) + 
		geom_bar(stat='identity',position ='dodge')+
		ylab(xlabel) + 
		xlab('Samples') + 
		scale_x_discrete(breaks = x$sample, label = str_replace(x$sample,'_(glucose|galactose)','')) + 
		theme(axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
		coord_flip() + 
		scale_y_log10() + 
		scale_fill_discrete(name = 'Treatment')

	return(p)
}



fn_list = list.files(star_alignment_dir, pattern = 'final.out', recursive=TRUE, full.names = TRUE)
fn_list = fn_list[!str_detect(fn_list,'set8/')]

stats = foreach(fn = fn_list, .combine = rbind)%do%{
	sample = basename(fn)
	sample = str_split_fixed(sample,'\\.Log',2)[,1]
	sample = str_replace(sample,'.R1.fastq.gz','')
	sample = tolower(sample)
	stats = read_STAR_mapping_stats(fn)
	stats$sample = sample
	return(stats)
}
stats[,condition := str_extract(sample,'(glucose|galactose)')]

input_reads = stats[field %in% c('Number of input reads')]
input_reads[, value := as.integer(value)]
summary(input_reads$value)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 24382915 45533884 52694050 52739197 60078076 92793601 
sum(as.numeric(input_reads$value))
# 2531481442
p1 = plot_STAR_mapping_stats(input_reads,'Number of input reads') + scale_fill_discrete(guide = 'none')

uniquely_mapped = stats[field %in% c('Uniquely mapped reads number')]
uniquely_mapped[, value := as.integer(value)]
summary(uniquely_mapped$value)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 18940303 41026853 46768544 47417802 55229122 84834862 
sum(as.numeric(uniquely_mapped$value))
p2 = plot_STAR_mapping_stats(uniquely_mapped,'Uniquely mapped reads')


p = plot_grid(p1,p2,nrow=1,rel_widths=c(1,1.4),labels = c('a','b'))
fig_fn = sprintf('%s/mapping_stats_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8,base_height=7)

p = plot_grid(p1,p2,nrow=1,rel_widths=c(1,1.4),labels = c('A','B'))
fig_fn = sprintf('%s/mapping_stats_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8,base_height=7)