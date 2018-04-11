library(data.table)
library(stringr)
library(gplots)

# Variables:
chr1_dosage_fn = '../processed_data/genotype_qc/detect_duplication/chr1_dosage.txt'
fig_dir = '../figures/genotype_qc/detect_duplication/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

# Functions:
read_dosage = function(fn){
	dosage = fread(fn)
	col_names = colnames(dosage)
	col_names = str_split_fixed(col_names, ']', 2)[,2]
	col_names = str_replace(col_names, ':DS', '')
	colnames(dosage) = col_names
	row_data = dosage[,1:4]
	dosage = dosage[,5:ncol(dosage)]
	return(list(row_data,dosage))
}

remove_multiallelic_sites = function(row_data,dosage){
	is_multiallelic = str_detect(row_data$ALT,',')
	row_data = row_data[!is_multiallelic]
	dosage = dosage[!is_multiallelic]
	return(list(row_data,dosage))
}

update_sample_names = function(old,old_to_new_fn = '../data/meta/dna2rna.txt'){
	old_to_new = fread(old_to_new_fn,colClasses = c('character','character'))
	setnames(old_to_new, c('old','new'))
	map = old_to_new$new
	names(map) = old_to_new$old
	return(map[old])
}

calculate_dosage_correlation = function(dosage){
	if (!is.matrix(dosage)){
		dosage = as.matrix(dosage)
	}
	mode(dosage) = 'numeric'
	cor = cor(dosage)
	return(cor)
}

plot_dosage_correlation = function(cor){
	heatmap.2(cor,Rowv = 'none', Colv = 'none', trace = 'none', dendrogram = 'none',margin=c(7,7))
}

# Main:
chr1_dosage = read_dosage(chr1_dosage_fn)
chr1_dosage = remove_multiallelic_sites(chr1_dosage[[1]],chr1_dosage[[2]])
row_data = chr1_dosage[[1]]
chr1_dosage = chr1_dosage[[2]]
colnames(chr1_dosage) = update_sample_names(colnames(chr1_dosage))
chr1_dosage = chr1_dosage[,order(colnames(chr1_dosage)),with=FALSE]
cor = calculate_dosage_correlation(chr1_dosage)

fig_fn = sprintf('%s/genotype_correlation.pdf',fig_dir)
pdf(fig_fn)
plot_dosage_correlation(cor)
dev.off()