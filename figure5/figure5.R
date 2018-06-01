library(data.table)
library(cowplot)
library(manhattan)
library(ggrepel)
library(stringr)
detach('package:locuscomparer',unload=TRUE)
devtools::install_github('boxiangliu/locuscomparer')
library(locuscomparer)

#-----------#
# Variables #
#-----------#
amd_eQTL_manhattan_rds = '../processed_data/finemap/manhattan/manhattan_amd/manhattan.rds'
amd_sQTL_manhattan_rds = '../processed_data/finemap/manhattan/manhattan_amd_sqtl/manhattan.rds'
myopia_eQTL_manhattan_rds = '../processed_data/finemap/manhattan/manhattan_myopia/manhattan.rds'
myopia_sQTL_manhattan_rds = '../processed_data/finemap/manhattan/manhattan_myopia_sqtl/manhattan.rds'
RDH5_eQTL_fn = '../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135437.5_RDH5.txt'
chr12_sQTL_fn = '../processed_data/sqtl/fastQTL/nominal/glucose/chr12.nominal.txt.gz'

fig_dir = '../figures/figure5/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
CUTOFF = 0.05
#-----------#
# Functions #
#-----------#
make_manhattan_data = function(eQTL_rds,sQTL_rds,cutoff = CUTOFF){
	eQTL = readRDS(eQTL_rds)
	eQTL_data = eQTL[[1]]
	eQTL_data$gene_id = NULL
	eQTL_data[,y:=abs(y)]

	sQTL = readRDS(sQTL_rds)
	sQTL_data = sQTL[[1]]
	sQTL_data$clu_name = NULL
	sQTL_data[,y:=abs(y)]

	eQTL_data = eQTL_data[condition == 'Glucose']
	sQTL_data = sQTL_data[condition == 'Glucose']

	eQTL_data$type = 'eQTL'
	sQTL_data$type = 'sQTL'

	manhattan_data = rbind(eQTL_data,sQTL_data)
	manhattan_data[,type:=factor(type,level=c('eQTL','sQTL'))]
	manhattan_data[,label:=ifelse(abs(y)>cutoff,gene_name,'')]

	return(manhattan_data)
}

plot_manhattan = function(data,label, cutoff = CUTOFF){
	dummy=data.table(type=c('eQTL','sQTL'),y=c(cutoff,cutoff))
	text_data = data.frame(cumulative_pos = 5e8, y = 0.4, color = 'black', type = factor('eQTL',levels = c('eQTL','sQTL')))
	p=manhattan(data,build='hg19')+
		facet_grid(type~.,scale='free_y')+
		scale_y_sqrt(breaks = c(0.05,seq(0.1,0.5,0.1)))+
		geom_text_repel(aes(label=label),color='black')+
		geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
		ylab(paste('Colocalization probability')) +
		theme_bw() + 
		theme(panel.grid = element_blank(), panel.border = element_rect(size = 1), strip.background = element_rect(color = 'black',fill = 'white', size = 1)) + 
		geom_text(data = text_data, label = label, size = 5)
	return(p)
}

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


read_amd_gwas = function(){
	gwas = fread('../data/gwas/Fritsche_2015_AdvancedAMD.txt')
	gwas[,Chrom:=paste0('chr',Chrom)]
	gwas = gwas[,list(rsid=Marker,chr=Chrom,pos=Pos,pval=GC.Pvalue)]
	return(gwas)
}

read_myopia_gwas = function(){
	gwas = fread('zcat ../data/gwas/23andme_myopia.prepared.txt.gz')
	gwas = gwas[,list(rsid,chr,pos=snp_pos,pval=pvalue)]
	return(gwas)
}

#------#
# Main #
#------#

# Manhattan:
amd_manhattan_data = make_manhattan_data(amd_eQTL_manhattan_rds,amd_sQTL_manhattan_rds)
amd_manhattan_plot = plot_manhattan(amd_manhattan_data,label='AMD') 

myopia_manhattan_data = make_manhattan_data(myopia_eQTL_manhattan_rds,myopia_sQTL_manhattan_rds)
myopia_manhattan_plot = plot_manhattan(myopia_manhattan_data, label = 'Myopia')

# Conservation:
conservation = plot_conservation()

# Locuscompare:
amd_gwas = read_amd_gwas()
myopia_gwas = read_myopia_gwas()
RDH5_eQTL = read_rasqual(RDH5_eQTL_fn)
RDH5_sQTL = read_fastQTL(chr12_sQTL_fn)[fid == 'chr12:56115278:56117670:clu_1202'][,c('chr','pos'):=extract_chr_pos(snp)]


RDH5_eQTL_into_amd = merge(RDH5_eQTL,amd_gwas[,list(chr,pos,rsid)],by=c('chr','pos'))
eQTL_amd_locuscompare = main(
	in_fn1 = amd_gwas[,list(rsid,pval)],
	in_fn2 = RDH5_eQTL_into_amd[,list(rsid,pval)],
	title1 = 'AMD',
	title2 = 'eQTL',
	vcf_fn = '/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
	combine = FALSE,
	legend = FALSE)
eQTL_amd_locuscompare = eQTL_amd_locuscompare$locuscompare + theme(axis.title = element_text(size=11))

RDH5_sQTL_into_amd = merge(RDH5_sQTL,amd_gwas[,list(chr,pos,rsid)],by=c('chr','pos'))
sQTL_amd_locuscompare = main(
	in_fn1 = amd_gwas[,list(rsid,pval)],
	in_fn2 = RDH5_sQTL_into_amd[,list(rsid,pval)],
	title1 = 'AMD',
	title2 = 'sQTL',
	vcf_fn = '/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
	combine = FALSE,
	legend = FALSE)
sQTL_amd_locuscompare = sQTL_amd_locuscompare$locuscompare + theme(axis.title = element_text(size=11))

RDH5_eQTL_into_myopia = merge(RDH5_eQTL,myopia_gwas[,list(chr,pos,rsid)],by=c('chr','pos'))
eQTL_myopia_locuscompare = main(
	in_fn1 = myopia_gwas[,list(rsid,pval)],
	in_fn2 = RDH5_eQTL_into_myopia[,list(rsid,pval)],
	title1 = 'Myopia',
	title2 = 'eQTL',
	vcf_fn = '/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
	combine = FALSE,
	legend = FALSE)
eQTL_myopia_locuscompare = eQTL_myopia_locuscompare$locuscompare + theme(axis.title = element_text(size=11))

RDH5_sQTL_into_myopia = merge(RDH5_sQTL,myopia_gwas[,list(chr,pos,rsid)],by=c('chr','pos'))
sQTL_myopia_locuscompare = main(
	in_fn1 = myopia_gwas[,list(rsid,pval)],
	in_fn2 = RDH5_sQTL_into_myopia[,list(rsid,pval)],
	title1 = 'Myopia',
	title2 = 'sQTL',
	vcf_fn = '/srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
	combine = FALSE,
	legend = FALSE)
sQTL_myopia_locuscompare = sQTL_myopia_locuscompare$locuscompare + theme(axis.title = element_text(size=11))
# save.image('figure5/figure5.rda')
load('figure5/figure5.rda')



# Make Figure 5:
top_right = plot_grid(eQTL_amd_locuscompare,sQTL_amd_locuscompare,nrow=2,align='v',labels=c('B','C'))
top = plot_grid(amd_manhattan_plot,top_right,nrow=1,align='h',labels = c('A',''),rel_widths = c(3,1))
bottom_right = plot_grid(eQTL_myopia_locuscompare,sQTL_myopia_locuscompare,nrow=2,align='v',labels=c('E','F'))
bottom = plot_grid(myopia_manhattan_plot,bottom_right,nrow=1,align = 'h',labels = c('C',''),rel_widths = c(3,1))
entire = plot_grid(top,bottom,nrow=2,align='v',labels = '')
fig_fn = sprintf('%s/figure5.pdf',fig_dir)
save_plot(fig_fn,entire,base_width=8,base_height=8)