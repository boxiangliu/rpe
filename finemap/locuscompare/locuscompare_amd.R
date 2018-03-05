library(data.table)
detach('package:locuscomparer',unload=TRUE)
install.packages('/srv/persistent/bliu2/locuscomparer/',source=TRUE,repos=NULL)
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

fig_dir = '../figures/finemap/locuscompare/locuscompare_amd/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)
glucose_list = c(
	'../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135437.5_RDH5.txt',
	'../processed_data/rasqual/output/glucose/joint/chr7/ENSG00000059378.8_PARP12.txt',
	'../processed_data/rasqual/output/glucose/joint/chr9/ENSG00000196363.5_WDR5.txt',
	'../processed_data/rasqual/output/glucose/joint/chr18/ENSG00000082397.11_EPB41L3.txt'
	)
galactose_list = c(
	'../processed_data/rasqual/output/galactose/joint/chr12/ENSG00000135437.5_RDH5.txt',
	'../processed_data/rasqual/output/galactose/joint/chr7/ENSG00000059378.8_PARP12.txt',
	'../processed_data/rasqual/output/galactose/joint/chr18/ENSG00000082397.11_EPB41L3.txt'
	)


gwas = fread('../data/gwas/Fritsche_2015_AdvancedAMD.txt')
gwas[,Chrom:=paste0('chr',Chrom)]

dir.create('../figures/finemap/locuscompare/locuscompare_amd/glucose/')
dir.create('../figures/finemap/locuscompare/locuscompare_amd/galactose/')
for (in_fn in glucose_list){
	print(in_fn)
	chr = str_extract(in_fn,'chr[0-9]+')
	fid = str_extract(in_fn,'ENSG.+(?=\\.txt$)')
	rasqual = read_rasqual(in_fn)
	rasqual = merge(rasqual,gwas[,list(chr=Chrom,pos=Pos,rsid=Marker)],by=c('chr','pos'))
	p = main(
		in_fn1 = rasqual[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid=Marker,pval=GC.Pvalue)],
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
	rasqual = merge(rasqual,gwas[,list(chr=Chrom,pos=Pos,rsid=Marker)],by=c('chr','pos'))
	p = main(
		in_fn1 = rasqual[,list(rsid,pval)],
		in_fn2 = gwas[,list(rsid=Marker,pval=GC.Pvalue)],
		vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
	)
	save_plot(
		filename = sprintf('%s/galactose/%s.pdf',fig_dir,fid),
		plot = p,
		base_height = 4,
		base_width = 8)
}
