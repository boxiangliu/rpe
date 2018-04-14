library(data.table)
library(foreach)
library(stringr)
library(cowplot)

verifyBamID_dir = '../processed_data/rnaseq_qc/verifyBamID/'
fig_dir = '../figures/rnaseq_qc/plot_verifyBamID/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_verifyBamID = function(fn){
	result = fread(fn)[,c(1,3,7,12)]
	sample = basename(fn)
	sample = str_split_fixed(sample,'_dedup',2)[,1]
	sample = tolower(sample)
	result$sample = sample
	return(result)
}

plot_verifyBamID = function(x){
	setorder(x,CHIPMIX)
	x[,sample := factor(sample, unique(sample))]
	ggplot(x, aes(sample, CHIPMIX)) + 
		geom_point() + 
		coord_flip() + 
		xlab('')+
		ylab('Proportion of contamination')
}

fn_list = list.files(verifyBamID_dir,pattern = 'selfSM', full.names = TRUE)

result = foreach(fn = fn_list,.combine='rbind')%do%{
	read_verifyBamID(fn)
}

p = plot_verifyBamID(result)
fig_fn = sprintf('%s/proportion_of_contamination.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8, base_height = 8)


fn_list = list.files(verifyBamID_dir,pattern = 'bestSM', full.names = TRUE)

result = foreach(fn = fn_list,.combine='rbind')%do%{
	read_verifyBamID(fn)
}

result[,list(`#SEQ_ID` == CHIP_ID)]
# all matched