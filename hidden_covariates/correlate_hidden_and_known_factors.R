library(data.table)
library(stringr)
library(openxlsx)
library(ggcorrplot)
library(cowplot)

hidden_covariates_dir = '../processed_data/hidden_covariates/sva/'
joint_hidden_covariates_fn = paste0(hidden_covariates_dir,'glugal_LibSizCorr_SVAFact7_MeanGeq10_ZeroLeq20.txt')
glucose_hidden_covariates_fn = paste0(hidden_covariates_dir,'glucose_LibSizCorr_SVAFact4_MeanGeq10_ZeroLeq20.txt')
galactose_hidden_covariates_fn = paste0(hidden_covariates_dir,'galactose_LibSizCorr_SVAFact5_MeanGeq10_ZeroLeq20.txt')
known_covariates_fn = '../processed_data/hidden_covariates/known_covariates/tables.xlsx'
fig_dir = '../figures/hidden_covariates/correlate_hidden_and_known_factors/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_hidden_covariates = function(fn){
	hidden_covariates = fread(fn)
	setnames(hidden_covariates,'rn','ID')
	hidden_covariates[,ID := str_replace(ID,'\\.','_')]
	colnames(hidden_covariates) = str_replace(colnames(hidden_covariates),'f','SV')
	return(hidden_covariates)
}

read_known_covariates = function(fn,batch_as_indicator = FALSE){
	x1 = read.xlsx(fn,sheet = 1, rows = 1:25)
	x2 = read.xlsx(fn,sheet = 2)
	setDT(x1);setDT(x2)
	x2[,ID := tolower(Sample)]
	x2[,Sample := str_split_fixed(ID,'_',2)[,1]]
	x = merge(x1,x2,by='Sample')
	x$Sample = NULL
	x[,Sex := as.integer(as.factor(Sex))]
	x[,Treatment := as.integer(as.factor(Treatment))]
	x[,Ancestry := NULL]
	if (batch_as_indicator){
		x[,Batch := as.character(Batch)]
		x[,Batch := factor(Batch,as.character(1:12))]
		batch_matrix = model.matrix(~Batch+0,data = x[,list(Batch)])
		x = cbind(x,data.table(batch_matrix))
		x$Batch = NULL
	} else {
		x[,Batch := as.integer(Batch)]
	}
	return(x)
}


plot_correlation = function(hidden,known,x_axis,y_axis){
	all_covariates = merge(hidden,known,by='ID')
	all_covariates_mat = as.matrix(all_covariates[,2:ncol(all_covariates)])
	mode(all_covariates_mat) = 'numeric'
	cor = cor(all_covariates_mat)
	cor = cor[x_axis,y_axis]
	ggcorrplot(cor, lab = TRUE)
}

joint_hidden_covariates = read_hidden_covariates(joint_hidden_covariates_fn)
known_covariates = read_known_covariates(known_covariates_fn)

x_axis = c('Sex',paste0('Genotype.PC',1:3),'Treatment','RIN','Batch')
y_axis = paste0('SV',1:7)
p1 = plot_correlation(joint_hidden_covariates,known_covariates,x_axis,y_axis)

glucose_hidden_covariates = read_hidden_covariates(glucose_hidden_covariates_fn)
x_axis = c('Sex',paste0('Genotype.PC',1:3),'RIN','Batch')
y_axis = paste0('SV',1:4)
p2 = plot_correlation(glucose_hidden_covariates,known_covariates,x_axis,y_axis)

galactose_hidden_covariates = read_hidden_covariates(galactose_hidden_covariates_fn)
x_axis = c('Sex',paste0('Genotype.PC',1:3),'RIN','Batch')
y_axis = paste0('SV',1:5)
p3 = plot_correlation(galactose_hidden_covariates,known_covariates,x_axis,y_axis)

blank = ggplot() + geom_blank()
row1 = plot_grid(blank,p1,blank,rel_widths = c(0.5,1,0.5), labels = c('','A',''), align = 'h', nrow = 1)
row2 = plot_grid(p2,p3,rel_widths = c(1,1), labels = c('B','C'), align = 'h', nrow = 1)
p = plot_grid(row1,row2,nrow=2,rel_heights = c(7,6))
fig_fn = sprintf('%s/factor_correlation.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8, base_height = 8)
