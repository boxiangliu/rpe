# Read RASQUAL result: 
tmp_file=tempfile()
system(sprintf("cat ../processed_data/rasqual/output/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
ras=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file)
ras[,fid:=str_split_fixed(fid,'_',2)[,1]]
ras[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
pdf('../figures/rasqual/pvalue_distribution.pdf')
hist(ras$pval,breaks=1000,main='RASQUAL: MAF > 0.05',xlab='p-value')
dev.off()