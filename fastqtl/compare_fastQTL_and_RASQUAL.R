library(data.table)
library(stringr)
library(cowplot)
library(dplyr)

# Read FastQTL result: 
fas=fread('zcat ../processed_data/fastqtl/eqtl/nominal_pass.glucose.txt.gz',col.names=c('fid','sid','dist','pval','beta','se'))

# Read RASQUAL result: 
tmp_file=tempfile()
system(sprintf("cat ../processed_data/rasqual/test/total/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
ras_total=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
system(sprintf("cat ../processed_data/rasqual/test/ase/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
ras_ase=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
system(sprintf("cat ../processed_data/rasqual/test/joint/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
ras_joint=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file)


# Take the ensembl ID as fid:
ras_total[,fid:=str_split_fixed(fid,'_',2)[,1]]
ras_ase[,fid:=str_split_fixed(fid,'_',2)[,1]]
ras_joint[,fid:=str_split_fixed(fid,'_',2)[,1]]


# Calculate p-value from chisq stat: 
ras_total[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
ras_ase[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
ras_joint[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


# Remove variants without RS id:
fas=fas[sid!='.']


# Merge FastQTL and RASQUAL result: 
merged=merge(fas%>%select(fid,sid,pval_fas=pval),ras_total%>%select(fid,sid,pval_ras_total=pval),by=c('fid','sid'))
merged=merge(merged,ras_ase%>%select(fid,sid,pval_ras_ase=pval),by=c('fid','sid'))
merged=merge(merged,ras_joint%>%select(fid,sid,pval_ras_joint=pval),by=c('fid','sid'))
merged[,c('logpval_fas','logpval_ras_total','logpval_ras_ase','logpval_ras_joint'):=list(-log10(pval_fas),-log10(pval_ras_total),-log10(pval_ras_ase),-log10(pval_ras_joint))]


# Make scatterplot:
pdf('../figures/fastqtl/compare_fastQTL_with_RASQUAL.pdf',height=16,width=16)
pairs(~logpval_fas+logpval_ras_total+logpval_ras_ase+logpval_ras_joint,merged)
dev.off()
png('../figures/fastqtl/compare_fastQTL_with_RASQUAL.png')
pairs(~logpval_fas+logpval_ras_total+logpval_ras_ase+logpval_ras_joint,merged)
dev.off()


# Find p-value threshold in total count model below which joint model would not be significant: 
merged[pval_ras_total>0.1&pval_ras_joint<0.05]
nrow(merged)
nrow(merged[pval_ras_total<0.05])
nrow(merged[pval_ras_joint<0.05])
nrow(merged[pval_ras_joint<0.01])

cutoff=seq(0.1,1,0.1)
num_sig=integer(length(cutoff))
for (i in seq_along(cutoff)){
	num_sig[i]=nrow(merged[pval_ras_total<=cutoff[i]&pval_ras_joint<=0.05])
	pct_sig=(num_sig/num_sig[length(num_sig)])*100
}
data=data.frame(num_sig=num_sig,pct_sig=pct_sig,cutoff=cutoff)
p=ggplot(data,aes(x=as.character(cutoff),y=num_sig,label=sprintf("%.2f%%",pct_sig)))+geom_bar(stat='identity')+geom_text(nudge_y=300)+xlab('Total-read model p-value cutoff')+ylab('eQTL in joint model (p-value < 0.05) ')
save_plot('../figures/fastqtl/num_sig_eqtl_vs_pval_cutoff.pdf',p,base_height=8,base_width=9)