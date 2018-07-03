library(data.table)
library(cowplot)
library(stringr)

in_fn='/srv/persistent/bliu2/rpe/processed_data/sex/preprocess_vcf/rpe.chr1_and_X.tsv'
fig_dir='../figures/sex/'
out_dir='../processed_data/sex/sex/'
meta_fn='../data/meta/dna2rna.txt'

if (!dir.exists(out_dir)){dir.create(out_dir)}


dosage=fread(in_fn,header=T)
het_b=data.table(dosage[CHROM!='X',5:ncol(dosage),with=F]==1)
het_X=data.table(dosage[CHROM=='X',5:ncol(dosage),with=F]==1)


het_mean_b=colMeans(het_b,na.rm=T)
het_mean_X=colMeans(het_X,na.rm=T)


stopifnot(all(names(het_mean_b)==names(het_mean_X)))
data=data.table(sample=names(het_mean_b),chr1=het_mean_b,chrX=het_mean_X)
data[,gender:=ifelse(chrX/chr1<0.3,'M','F')]

meta = fread(meta_fn,header=TRUE,colClasses = c('character','character','integer'))
data = merge(data,meta[,list(DNA,NAME)],by.x='sample',by.y='DNA')

setorder(data,NAME)
data[,sample:=str_replace(sample,'\\.','-')]
data[,sample:=factor(sample,levels=rev(sample))]
data[,NAME:=factor(NAME,levels=rev(NAME))]
fwrite(data[,.(sample,gender)],sprintf('%s/gender.tsv',out_dir),sep='\t')


to_plot=melt(data,id.var=c('NAME','gender'),measure.vars=c('chr1','chrX'),variable.name='type',value.name='het_rate')

p1=ggplot(to_plot,aes(het_rate,NAME,color=type,label=gender))+
	geom_point(size=5)+
	geom_text(x=0.3,color='black')+
	background_grid(major = "xy", minor = "none")+
	xlab('Proportion Heterozygous')+
	ylab('')+
	scale_color_discrete(name = 'Chromosome')+
	xlim(0,0.3)
save_plot(sprintf('%s/proportion_heterozygous.pdf',fig_dir),p1,base_width=8,base_height=8)
