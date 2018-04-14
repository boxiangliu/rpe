library(data.table)
library(stringr)
library(openxlsx)
library(ggcorrplot)

hidden_covariates_fn = '../processed_data/hidden_covariates/sva/glugal_LibSizCorr_SVAFact7_MeanGeq10_ZeroLeq20.txt'
known_covariates_fn = '../processed_data/hidden_covariates/known_covariates/tables.xlsx'
read_hidden_covariates = function(fn){
	hidden_covariates = fread(fn)
	setnames(hidden_covariates,'rn','ID')
	hidden_covariates[,ID := str_replace(ID,'\\.','_')]
	return(hidden_covariates)
}

read_known_covariates = function(fn){
	x1 = read.xlsx(fn,sheet = 1, rows = 1:25)
	x2 = read.xlsx(fn,sheet = 2)
	setDT(x1);setDT(x2)
	x2[,ID := tolower(Sample)]
	x2[,Sample := str_split_fixed(ID,'_',2)[,1]]
	x = merge(x1,x2,by='Sample')
	return(x)
}
hidden_covariates = read_hidden_covariates(hidden_covariates_fn)
known_covariates = read_known_covariates(known_covariates_fn)

all_covariates = merge(hidden_covariates,known_covariates,by='ID')
all_covariates[,Sex := as.integer(as.factor(Sex))]
all_covariates[,Treatment := as.integer(as.factor(Treatment))]
all_covariates[,Ancestry := as.integer(as.factor(Ancestry))]
all_covariates_mat = as.matrix(all_covariates[,c(2:8,10:14)])
mode(all_covariates_mat) = 'numeric'
cor = cor(all_covariates_mat)
ggcorrplot(cor, lab = TRUE, type = "lower")