library(data.table)
library(locuscomparer)
library(stringr)

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

fig_dir = '../figures/finemap/locuscompare/locuscompare_amd_sqtl/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

glucose_list = c(
	RDH5 = 'chr12:56115278:56117670:clu_1202'
	)
galactose_list = c(
	)


gwas = fread('../data/gwas/Fritsche_2015_AdvancedAMD.txt')[,Chrom:=paste0('chr',Chrom)]
glucose_sqtl = read_fastQTL('../processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz')
setkey(glucose_sqtl,fid)

for (i in seq_along(glucose_list)){
	gene_name = names(glucose_list)[i]
	clu_name = glucose_list[i]
	print(gene_name)
	fastQTL = glucose_sqtl[fid == clu_name,][,c('chr','pos'):=extract_chr_pos(snp)]
	fastQTL = merge(fastQTL,gwas[,list(chr=Chrom,pos=Pos,rsid=Marker)],by=c('chr','pos'))
	chr = unique(fastQTL$chr)
	p = main(
		in_fn1 = fastQTL[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid=Marker,pval=GC.Pvalue)],
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
	fastQTL = merge(fastQTL,gwas[,list(chr=Chrom,pos=Pos,rsid=Marker)],by=c('chr','pos'))
	chr = unique(fastQTL$chr)
	p = main(
		in_fn1 = fastQTL[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid=Marker,pval=GC.Pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/galactose/%s_%s.pdf',fig_dir,gene_name,clu_name),
		plot = p,
		base_height = 4,
		base_width = 8)
}

