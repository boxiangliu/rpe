# Read in, format, and write out cell-based assay phenotypes for 24 hfRPE cell lines

library(xlsx)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)

phenotypes = c("id","trans_epithelial_resistance", "pigmentation", "nonmitochontrial_respiration","proton_leak","basal_respiration", "atp_production","maximal_respiration", "spare_reserve_capacity", "nonglycolytic_acidification","glycolysis","glycolytic_capacity","glycolytic_reserve","galactose_bound_os","galactose_total_os","galactose_ingested_os","glucose_bound_os","glucose_total_os", "glucose_ingested_os","regular_bound_os","regular_total_os", "regular_ingested_os")
excel_fn = '../processed_data/phenotype_association/format_phenotype/hfRPE_phenotypes.xlsx'
total = 24
out_dir = '../processed_data/phenotype_association/format_phenotype/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_excel_sheet = function(excel_fn, sheetIndex){
	data = read.xlsx(excel_fn, sheetIndex = i)
	data = data[c(-1,-2),]
	colnames(data) = phenotypes

	id = as.character(data$id[1])
	data$id = NULL

	melt_data = melt(data, id.vars = character(0))
	melt_data = data.table(melt_data)

	melt_data = melt_data[!is.na(value),]
	melt_data$sample = id

	return(melt_data)
}

data = foreach(i = 1:total, .combine = 'rbind')%do%{
	print(i)
	read_excel_sheet(excel_fn,i)
}

data[sample == '020206', sample := '020206-2']

setcolorder(data,c('sample','variable','value'))
out_fn = sprintf('%s/phenotypes.txt',out_dir)
fwrite(data,out_fn)