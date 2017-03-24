library(data.table)
library(stringr)
library(cowplot)

# Read FastQTL result: 
fas=fread('zcat ../processed_data/fastqtl/eqtl/nominal_pass.glucose.txt.gz',col.names=c('fid','sid','dist','pval','beta','se'))

# Read RASQUAL result: 
tmp_file=tempfile()
system(sprintf("cat ../processed_data/rasqual/output/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
ras=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file)


# Take the ensembl ID as fid:
ras[,fid:=str_split_fixed(fid,'_',2)[,1]]


# Calculate p-value from chisq stat: 
ras[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


# Merge FastQTL and RASQUAL result: 
merged=merge(fas,ras,by=c('fid','sid'))
merged_pval=merged%>%select(fid,sid,pval_fas=pval.x,pval_ras=pval.y)


# Make scatterplot:
p1=ggplot(merged_pval[1:100000],aes(pval_fas,pval_ras))+geom_point()+xlab('FastQTL p-value')+ylab('RASQUAL p-value')
p2=ggplot(merged_pval,aes(-log10(pval_fas),-log10(pval_ras)))+geom_point()+xlab('FastQTL -log10(p-value)')+ylab('RASQUAL -log10(p-value)')
pdf('../figures/fastqtl/compare_fastQTL_with_RASQUAL.pdf')
print(p1);print(p2);
dev.off()