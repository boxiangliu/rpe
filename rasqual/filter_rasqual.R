library(data.table)
library(cowplot)
library(R.utils)

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'

phi_ub = 0.75; phi_lb = 0.25
delta_ub = 0.1 
r2_fsnp_lb = r2_rsnp_lb = 0.9

tmp_file=tempfile()
system(sprintf("cat %s/*/*.txt > %s",glucose_dir,tmp_file))
glucose=fread(tmp_file,
	col.names=c('fid','sid','chr','pos','ref','alt','af',
		'hwe_chisq','impute_score','qval','chisq','pi',
		'delta','phi','overdispersion','snp_id',
		'num_fsnp','num_tsnp','n_iteration_null','n_iteration_alt',
		'ties','loglik','converge','r2_fsnp','r2_rsnp'))
file.remove(tmp_file)

glucose_filt = glucose[phi < phi_ub & phi > phi_lb & 
	delta < delta_ub & r2_fsnp > r2_fsnp_lb & 
	r2_rsnp > r2_rsnp_lb, ]

out_fn = sprintf('%s/all_associations_filt.txt', glucose_dir)
fwrite(glucose_filt, out_fn)
gzip(out_fn)


tmp_file=tempfile()
system(sprintf("cat %s/*/*.txt > %s",galactose_dir,tmp_file))
galactose=fread(tmp_file,
	col.names=c('fid','sid','chr','pos','ref','alt','af',
		'hwe_chisq','impute_score','qval','chisq','pi',
		'delta','phi','overdispersion','snp_id',
		'num_fsnp','num_tsnp','n_iteration_null','n_iteration_alt',
		'ties','loglik','converge','r2_fsnp','r2_rsnp'))
file.remove(tmp_file)

galactose_filt = galactose[phi < phi_ub & phi > phi_lb & 
	delta < delta_ub & r2_fsnp > r2_fsnp_lb & 
	r2_rsnp > r2_rsnp_lb, ]

out_fn = sprintf('%s/all_associations_filt.txt', galactose_dir)
fwrite(galactose_filt, out_fn)
gzip(out_fn)
