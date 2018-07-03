library(data.table)
library(locuscomparer)
library(stringr)
library(cowplot)
library(ggrepel)

read_rasqual = function(in_fn){
	x = fread(
		input = in_fn,
		select = c(1,2,3,4,11,12),
		col.names = c('fid','snp','chr','pos','chisq','pi'))
	x[,gene_name := str_split_fixed(fid,'_',2)[,2]]
	x[,gene_id := str_split_fixed(fid,'_',2)[,1]]
	x[,pval := pchisq(q=chisq,df=1,lower.tail=FALSE)]
	return(x)
}

fig_dir = '../figures/finemap/locuscompare/locuscompare_myopia/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

glucose_list = c(
	'../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135437.5_RDH5.txt',
	'../processed_data/rasqual/output/glucose/joint/chr8/ENSG00000120885.15_CLU.txt',
	'../processed_data/rasqual/output/glucose/joint/chr15/ENSG00000138613.9_APH1B.txt',
	# '../processed_data/rasqual/output/glucose/joint/chr19/ENSG00000248098.6_BCKDHA.txt',
	'../processed_data/rasqual/output/glucose/joint/chr19/ENSG00000104805.11_NUCB1.txt',
	'../processed_data/rasqual/output/glucose/joint/chr14/ENSG00000187097.8_ENTPD5.txt',
	'../processed_data/rasqual/output/glucose/joint/chr2/ENSG00000240344.4_PPIL3.txt',
	'../processed_data/rasqual/output/glucose/joint/chr21/ENSG00000157557.7_ETS2.txt'
	)

galactose_list = c(
	'../processed_data/rasqual/output/galactose/joint/chr12/ENSG00000135437.5_RDH5.txt',
	'../processed_data/rasqual/output/galactose/joint/chr8/ENSG00000120885.15_CLU.txt',
	'../processed_data/rasqual/output/galactose/joint/chr12/ENSG00000172572.6_PDE3A.txt'
	)


gwas = fread('zcat ../data/gwas/23andme_myopia.prepared.txt.gz')

for (in_fn in glucose_list){
	print(in_fn)
	chr = str_extract(in_fn,'chr[0-9]+')
	fid = str_extract(in_fn,'ENSG.+(?=\\.txt$)')
	rasqual = read_rasqual(in_fn)
	rasqual = merge(rasqual,gwas[,list(chr,pos=snp_pos,rsid)],by=c('chr','pos'))
	p = main(
		in_fn1 = rasqual[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid,pval=pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/glucose/%s.pdf',fig_dir,fid),
		plot = p,
		base_height = 4,
		base_width = 8)
}

for (in_fn in galactose_list){
	print(in_fn)
	chr = str_extract(in_fn,'chr[0-9]+')
	fid = str_extract(in_fn,'ENSG.+(?=\\.txt$)')
	rasqual = read_rasqual(in_fn)
	rasqual = merge(rasqual,gwas[,list(chr,pos=snp_pos,rsid)],by=c('chr','pos'))
	p = main(
		in_fn1 = rasqual[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid,pval=pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/galactose/%s.pdf',fig_dir,fid),
		plot = p,
		base_height = 4,
		base_width = 8)
}
