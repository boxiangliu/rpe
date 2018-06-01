library(data.table)
library(foreach)

orig_dir = '../data/genotype/orig/qc/'
orig_fn_list = list.files(orig_dir,pattern='stats',full.names=TRUE)
filt_dir = '../data/genotype/filt/qc/'
filt_fn_list = list.files(filt_dir,pattern='stats',full.names=TRUE)

count_variant = function(fn_list){
	count = foreach(fn = fn_list,.combine='c')%dopar%{
		count = fread(sprintf('grep "number of records" %s',fn))
		count = count$V4
		return(count)
	}
	return(sum(count))
}


orig_count = count_variant(orig_fn_list) 
orig_count # 30761499
filt_count = count_variant(filt_fn_list)
filt_count # 12955723