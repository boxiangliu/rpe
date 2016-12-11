library(data.table)
library(dtplyr)
library(dplyr)
library(cowplot)
library(gap)
library(stringr)
source('gwas_atacseq_overlap/utils.R')

## Functions: 
define_peak_region=function(x,col=c('chr','start','end'),size=100){
	setnames(x,col,c('chr','start','end'))
	x[,peak_pos:=start+peak]
	x[,start:=peak_pos-size]
	x[,end:=peak_pos+size]
	x[,peak_pos:=NULL]
	setnames(x,c('chr','start','end'),col)
	return(x)
}


## Variables: 
HRPEpiC_file='/srv/persistent/bliu2/rpe/data/atacseq/HRPEpiC/GSM736630_hg19_wgEncodeUwDnaseHrpePkRep1V2.narrowPeak.gz'
gwas_file='../data/gwas/Fritsche_2015_AdvancedAMD.txt'
rpe_file='../data/atacseq/rpe/out/peak/idr/pseudo_reps/rep1/rpe_rep1-pr.IDR0.1.filt.narrowPeak.gz'
Roadmap_dir='/mnt/data/epigenomeRoadmap/peaks/consolidated/narrowPeak/'

# read GWAS data:
gwas=fread(gwas_file)


# format GWAS data:
gwas=gwas%>%select(chrom=Chrom,pos=Pos,rsid=Marker,pval=GC.Pvalue)
gwas=gwas%>%mutate(chromStart=pos,chromEnd=pos,chrom=paste0('chr',chrom))
gwas=gwas%>%select(chrom,chromStart,chromEnd,rsid,pval)



# read ATACseq peaks: 
atac_rpe=readAtac(paste0('zcat ',rpe_file))
atac_rpe=define_peak_region(atac_rpe,col=c('chrom','chromStart','chromEnd'),size=75)
atacseq=list(atac_rpe)


# Read Encode HRPEpiC sample: 
atacseq[[length(atacseq)+1]]=readAtac(sprintf('zcat %s',HRPEpiC_file))


# read roadmap atacseq samples:
DNase_samples=list.files(Roadmap_dir,pattern='*DNase.macs2.narrowPeak*',full.names=T)
for (i in 1:length(DNase_samples)){
	atacseq[[length(atacseq)+1]]=define_peak_region(readAtac(sprintf('zcat %s',DNase_samples[i])),col=c('chrom','chromStart','chromEnd'),size=75)
}



# assign tissue names:
tissues=c('rpe','HRPEpiC',str_match(DNase_samples,'E[0-9]{3}'))
names(atacseq)=tissues

# overlap summary: 
# thd=c(1e-3,1e-4,1e-5,1e-6,1e-7)
thd=c(1e-3)
baseline=nrow(atacseq[['rpe']])
overlap_summary=data.frame()
for (this_thd in thd){
	# thresholding gwas by pvalue: 
	gwas_thd=gwas[pval<this_thd,]

	for (this_tissue in names(atacseq)[1:2]){
	# for (this_tissue in names(atacseq)){
		# get atacseq data:
		atac=atacseq[[this_tissue]]


		# calculate overlaps: 
		overlaps=overlap(atac,gwas_thd)


		# Number of GWAS SNPs in each peak:
		assign(paste0('overlaps_peakSummary_',this_tissue),peakSummary(overlaps))


		# Number of unique peaks overlapping GWAS SNPs:
		n_unique_peaks=nrow(get(paste0('overlaps_peakSummary_',this_tissue)))


		# append to overlap summary: 
		tmp_df=data.frame(ratio=nrow(atac)/baseline,fraction=n_unique_peaks/nrow(gwas_thd),tissue=this_tissue,threshold=this_thd)
		overlap_summary=rbind(overlap_summary,tmp_df)
	}
}


# calculating fractions:
overlap_summary=as.data.table(overlap_summary)
# overlap_summary[,fraction_uniq:=n_unique_peaks/n_peaks]
# overlap_summary[,fraction:=n_overlaps/n_peaks]
overlap_summary[,fraction_adj:=fraction/ratio]
overlap_summary

# calculating lfc:
overlap_summary[,rank:=rank(fraction_uniq),by='threshold']
overlap_summary=overlap_summary%>%arrange(threshold,desc(rank))
lfc=log10(overlap_summary[rank==max(rank),fraction_uniq]/overlap_summary[rank==max(rank)-1,fraction_uniq])
annotation=data.frame(overlap_summary[rank==max(rank),.(threshold,fraction_uniq,tissue)],lfc=prettyNum(lfc,digits=3))


# plotting: 
p=ggplot(overlap_summary,aes(x=threshold,y=fraction_uniq,color=tissue,size=ifelse(tissue=='HRPEpiC',2,1)))+geom_point()+scale_x_log10(breaks=10^(-seq(3,7)))+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+xlab('Threshold')+ylab('Fraction of unique overlaps')
p=ggplot(overlap_summary,aes(x=threshold,y=fraction_adj*100,color=tissue,size=ifelse(tissue=='HRPEpiC',2,1)))+geom_point()+scale_x_log10(breaks=10^(-seq(3,7)))+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+xlab('Threshold')+ylab('Fraction of unique overlaps')
p
save_plot('../figures/gwas_atacseq_overlap/gwas_atacseq_overlap.unique.pdf',p,base_height=6,base_width=6)

p1=ggplot(overlap_summary[threshold==1e-7,],aes(x=reorder(tissue,-fraction_uniq,FUN=mean),y=fraction_uniq,color=tissue))+geom_point()+scale_y_log10(breaks=10^(-seq(2,4,by=0.3)))+theme_bw()+xlab('Tissue')+ylab('Fraction of unique overlaps')+theme(axis.text.x=element_text(angle=45,hjust=1))
save_plot('../figures/gwas_atacseq_overlap/gwas_atacseq_overlap.unique.1e-7.pdf',p1,base_height=6,base_width=8)

