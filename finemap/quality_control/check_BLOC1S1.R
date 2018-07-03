library(data.table)
library(stringr)
library(cowplot)
library(locuscomparer)
library(ggrepel)

bloc1s1_fn = '../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135441.3_BLOC1S1.txt'
fig_dir = '../figures/finemap/quality_control/check_BLOC1S1/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

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

bloc1s1 = read_rasqual(bloc1s1_fn)
gwas = fread('../data/gwas/Fritsche_2015_AdvancedAMD.txt')
gwas[,Chrom:=paste0('chr',Chrom)]
bloc1s1 = merge(bloc1s1,gwas[,list(chr=Chrom,pos=Pos,rsid=Marker)],by=c('chr','pos'))
chr = 'chr12'
p = main(
	in_fn1 = bloc1s1[,list(rsid,pval)],
	in_fn2 = gwas[,list(rsid=Marker,pval=GC.Pvalue)],
	vcf_fn = sprintf('/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',chr)
)
fig_fn = sprintf('%s/BLOC1S1.pdf',fig_dir)
save_plot(fig_fn,p,base_width=8,base_height=4)
