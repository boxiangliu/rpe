library(data.table)
library(doMC)
library(foreach)
registerDoMC(11)
library(stringr)
library(cowplot)
library(ggrepel)
library(gtools)

rpe_dir='../processed_data/genotype_pc/preprocess_vcf/rpe'
okg_dir='../processed_data/genotype_pc/preprocess_vcf/1kg'
out_dir='../processed_data/genotype_pc/genotype_pc/'
fig_dir='../figures/genotype_pc/genotype_pc/'
meta_fn='../data/meta/dna2rna.txt'
if(!dir.exists(out_dir)) dir.create(out_dir)
if(!dir.exists(fig_dir)) dir.create(fig_dir,recursive=T)

plot_screeplot = function(pc,top = 10){
	variance = pc$sdev^2
	total_variance = sum(variance)
	proportion = variance/total_variance
	data = data.table(pc = seq_along(proportion), proportion)
	ggplot(data[1:top], aes(pc, proportion)) + 
		geom_point() + 
		geom_line() + 
		xlab('PC') + 
		ylab('Proportion of variance') + 
		scale_x_continuous(breaks = 1:top, labels = 1:top)
}

# Read RPE genotype dosage:
in_fn=list.files(rpe_dir,pattern='chr')
in_fn=mixedsort(in_fn)
dosage=foreach(f=in_fn,.combine='rbind')%dopar%{fread(sprintf("%s/%s",rpe_dir,f))}


# Some reformatting for RPE: 
new_names=str_replace(str_split_fixed(names(dosage),']',2)[,2],':DS','')
setnames(dosage,new_names)
setDF(dosage)
rownames(dosage)=dosage$ID
dosage$ID=NULL


# PCA:
dosage_t=t(dosage)
pc=prcomp(dosage_t,center=TRUE,scale=FALSE)


# Save PCA results:
out=as.data.frame(pc$x)
out$sample=rownames(out)
setcolorder(out,c(ncol(out),1:(ncol(out)-1)))
fwrite(out,sprintf('%s/pc.tsv',out_dir),sep='\t')


# Plot the first few PCs: 
to_plot=as.data.frame(pc$x[,1:10])
to_plot$sample=rownames(to_plot)
pdf(sprintf("%s/pc.pdf",fig_dir))
plot(pc)
ggplot(to_plot,aes(PC1,PC2,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC1,PC3,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC2,PC3,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
dev.off()


# Read 1kg genotype dosage: 
in_fn=list.files(okg_dir,pattern='chr')
in_fn=mixedsort(in_fn)
dosage_1kg_chr1=foreach(f=in_fn[1],.combine='rbind')%dopar%{fread(sprintf("%s/%s",okg_dir,f))}


# Some reformatting for RPE: 
new_names=str_replace(str_split_fixed(names(dosage_1kg_chr1),']',2)[,2],':GT','')
setnames(dosage_1kg_chr1,new_names)


# Merge RPE and 1kg:
dosage$ID=rownames(dosage)
setDT(dosage)
dosage_merged=merge(dosage,dosage_1kg_chr1,by='ID',sort=FALSE)
setDF(dosage_merged)
rownames(dosage_merged)=dosage_merged$ID
dosage_merged$ID=NULL


# Perform PCA:
dosage_merged_t=t(dosage_merged)
pc=prcomp(dosage_merged_t,center=TRUE,scale=FALSE)


# Save PCA results:
out=as.data.frame(pc$x)
out$sample=rownames(out)
setcolorder(out,c(ncol(out),1:(ncol(out)-1)))
fwrite(out,sprintf('%s/pc_with_1kg.tsv',out_dir),sep='\t')


# Read in 1kg panel file:
panel=fread('../processed_data/genotype_pc/1kg_samples/4_sample_each_pop.panel')


# Plot the first few PCs: 
to_plot=as.data.frame(pc$x[,1:10])
to_plot$sample=rownames(to_plot)
to_plot=merge(to_plot,panel,by='sample',all=T,sort=FALSE)
setDT(to_plot)
to_plot[,super_pop:=ifelse(is.na(super_pop),'unknown',super_pop)]
to_plot[,label:=ifelse(super_pop=='unknown',sample,'')]


pdf(sprintf("%s/pc_with_1kg.pdf",fig_dir))
plot(pc)
ggplot(to_plot,aes(PC1,PC2,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC1,PC3,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC2,PC3,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
dev.off()


# Remove duplicate sample 021011:
dosage$`021011`=NULL


# PCA:
dosage_t=t(dosage)
pc=prcomp(dosage_t,center=TRUE,scale=FALSE)


# Save PCA results:
out=as.data.frame(pc$x)
out$sample=rownames(out)
setcolorder(out,c(ncol(out),1:(ncol(out)-1)))
fwrite(out,sprintf('%s/pc_nodup.tsv',out_dir),sep='\t')


# Plot the first few PCs: 
to_plot=as.data.frame(pc$x[,1:10])
to_plot$sample=rownames(to_plot)
pdf(sprintf("%s/pc_nodup.pdf",fig_dir))
plot(pc)
ggplot(to_plot,aes(PC1,PC2,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC1,PC3,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
ggplot(to_plot,aes(PC2,PC3,label=sample))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)
dev.off()


# Merge RPE and 1kg:
dosage_merged=merge(dosage,dosage_1kg_chr1,by='ID',sort=FALSE)
setDF(dosage_merged)
rownames(dosage_merged)=dosage_merged$ID
dosage_merged$ID=NULL


# Perform PCA:
dosage_merged_t=t(dosage_merged)
pc=prcomp(dosage_merged_t,center=TRUE,scale=FALSE)


# Save PCA results:
saveRDS(pc,sprintf('%s/pc_with_1kg_nodup.rds',out_dir))
Sys.chmod(sprintf('%s/pc_with_1kg_nodup.rds',out_dir), "444", use_umask = FALSE)
pc = readRDS(sprintf('%s/pc_with_1kg_nodup.rds',out_dir))
out=as.data.frame(pc$x)
out$sample=rownames(out)
setcolorder(out,c(ncol(out),1:(ncol(out)-1)))
fwrite(out,sprintf('%s/pc_with_1kg_nodup.tsv',out_dir),sep='\t')


# Read in 1kg panel file:
panel=fread('../processed_data/genotype_pc/1kg_samples/4_sample_each_pop.panel')


# Plot the first few PCs: 
to_plot=as.data.frame(pc$x[,1:10])
to_plot$sample=rownames(to_plot)

meta = fread(meta_fn,header=TRUE,colClasses = c('character','character','integer'))
to_plot$sample[1:23] = meta$NAME[match(to_plot$sample[1:23],meta$RNA)]

to_plot=merge(to_plot,panel,by='sample',all=T,sort=FALSE)
setDT(to_plot)
to_plot[,super_pop:=ifelse(is.na(super_pop),'unknown',super_pop)]
to_plot[,label:=ifelse(super_pop=='unknown',sample,'')]

pdf(sprintf("%s/screeplot_with_1kg_nodup.pdf",fig_dir))
plot_screeplot(pc)
dev.off()

pdf(sprintf("%s/pc_with_1kg_nodup.pdf",fig_dir))
plot(pc)
ggplot(to_plot,aes(PC1,PC2,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)+scale_color_discrete(name = 'Ancestry', breaks = c('AFR','AMR','EAS','EUR','SAS','unknown'), labels =c('African','Native American','East Asian','European','South Asian','Unknown (RPE)'))
ggplot(to_plot,aes(PC1,PC3,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)+scale_color_discrete(name = 'Ancestry', breaks = c('AFR','AMR','EAS','EUR','SAS','unknown'), labels =c('African','Native American','East Asian','European','South Asian','Unknown (RPE)'))
ggplot(to_plot,aes(PC2,PC3,color=super_pop,label=label))+geom_point(size=5,alpha=0.5)+geom_text_repel(force=3)+scale_color_discrete(name = 'Ancestry', breaks = c('AFR','AMR','EAS','EUR','SAS','unknown'), labels =c('African','Native American','East Asian','European','South Asian','Unknown (RPE)'))
dev.off()