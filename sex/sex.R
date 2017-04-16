library(data.table)
library(cowplot)

in_fn='/srv/persistent/bliu2/rpe/processed_data/sex/preprocess_vcf/rpe.chr1_and_X.tsv'
fig_dir='../figures/sex/'
out_dir='../processed_data/sex/sex/'
if (!dir.exists(out_dir)){dir.create(out_dir)}


dosage=fread(in_fn,header=T)
het_b=data.table(dosage[CHROM!='X',5:ncol(dosage),with=F]==1)
het_X=data.table(dosage[CHROM=='X',5:ncol(dosage),with=F]==1)


het_mean_b=colMeans(het_b,na.rm=T)
het_mean_X=colMeans(het_X,na.rm=T)


stopifnot(all(names(het_mean_b)==names(het_mean_X)))
data=data.table(sample=names(het_mean_b),autosome=het_mean_b,chrX=het_mean_X)
data[,gender:=ifelse(chrX/autosome<0.3,'M','F')]
fwrite(data[,.(sample,gender)],sprintf('%s/gender.tsv',out_dir),sep='\t')


to_plot=melt(data,id.var=c('sample','gender'),measure.vars=c('autosome','chrX'),variable.name='type',value.name='het_rate')
p1=ggplot(to_plot,aes(het_rate,sample,color=type,label=gender))+geom_point(size=5)+geom_text(x=0.3,color='black')+background_grid(major = "xy", minor = "none")+xlab('Proportion Heterozygous')+xlim(0,0.3)
save_plot(sprintf('%s/proportion_heterzygous.pdf',fig_dir),p1,base_width=8,base_height=8)


