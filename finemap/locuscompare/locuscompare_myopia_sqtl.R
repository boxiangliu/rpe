library(data.table)
library(locuscomparer)
library(stringr)
library(cowplot)
library(ggrepel)

read_fastQTL = function(in_fn){
	x = fread(
		input = sprintf('gunzip -c %s',in_fn),
		select = c(1,2,4),
		col.names = c('fid','snp','pval'))
	return(x)
}

extract_chr_pos = function(x){
	split_x = str_split_fixed(x,'_',5)
	chr = split_x[,1]
	pos = as.integer(split_x[,2])
	if (!str_detect(chr[1],'chr')){
		chr = paste0('chr',chr)
	}
	return(list(chr,pos))
}

fig_dir = '../figures/finemap/locuscompare/locuscompare_myopia_sqtl/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

glucose_list = c(
	RDH5 = 'chr12:56115278:56117670:clu_1202',
	# TCP11L2 = 'chr12:106706212:106708136:clu_1434',
	# YTHDC1 = 'chr4:69197915:69199026:clu_11700',
	# RNF13 = 'chr3:149532118:149563798:clu_12913',
	TSPAN10 = 'chr17:79612655:79614962:clu_4491'
	)
galactose_list = c(
	RDH5 = 'chr12:56115278:56117670:clu_1247',
	# RWDD2A = 'chr6:83904371:83904914:clu_10829',
	TAGLN = 'chr11:117070545:117073718:clu_2714',
	ANKRD9 = 'chr14:102974886:102975866:clu_6681',
	SCAMP3 = 'chr1:155230450:155231877:clu_16149',
	FBXO7 = 'chr22:32883815:32889092:clu_8550',
	PPIL3 = 'chr2:201741760:201742220:clu_14556'
	)


gwas = fread('zcat ../data/gwas/23andme_myopia.prepared.txt.gz')
glucose_sqtl = read_fastQTL('../processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz')
setkey(glucose_sqtl,fid)

for (i in seq_along(glucose_list)){
	gene_name = names(glucose_list)[i]
	clu_name = glucose_list[i]
	print(gene_name)
	fastQTL = glucose_sqtl[fid == clu_name,][,c('chr','pos'):=extract_chr_pos(snp)]
	fastQTL = merge(fastQTL,gwas[,list(chr,pos=snp_pos,rsid)],by=c('chr','pos'))
	chr = unique(fastQTL$chr)
	p = main(
		in_fn1 = fastQTL[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid,pval=pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/glucose/%s_%s.pdf',fig_dir,gene_name,clu_name),
		plot = p,
		base_height = 4,
		base_width = 8)
}

galactose_sqtl = read_fastQTL('../processed_data/sqtl/fastQTL/nominal/galactose/all.nominal.txt.gz')
setkey(galactose_sqtl,fid)
for (i in seq_along(galactose_list)){
	gene_name = names(galactose_list)[i]
	clu_name = galactose_list[i]
	print(gene_name)
	fastQTL = galactose_sqtl[fid == clu_name,][,c('chr','pos'):=extract_chr_pos(snp)]
	fastQTL = merge(fastQTL,gwas[,list(chr,pos=snp_pos,rsid)],by=c('chr','pos'))
	chr = unique(fastQTL$chr)
	p = main(
		in_fn1 = fastQTL[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid,pval=pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/galactose/%s_%s.pdf',fig_dir,gene_name,clu_name),
		plot = p,
		base_height = 4,
		base_width = 8)
}