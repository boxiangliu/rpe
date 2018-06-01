# Count the number of sQTLs:

library(data.table)
library(foreach)
library(stringr)

args=commandArgs(T)
in_fn=args[1]
out_dir=args[2]
out_prefix=args[3]
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

# Functions:
main=function(in_fn,out_dir,out_prefix,fdr=c(0.001,0.01,0.05)){
	in_fn = '../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz'
	
	sqtl = fread(sprintf('zcat %s',in_fn),select=c(1,6,7,16),col.names=c('intron','snp','dist','pval'))
	sqtl[, cluster := str_split_fixed(intron,':',4)[,4]]
	sqtl[,fwer:=p.adjust(pval,'bonferroni'),by='cluster']
	set.seed(42)
	sqtl[,fwer_rank:=rank(fwer,ties.method='random'),by='cluster']
	cluster_pval = sqtl[fwer_rank==1,list(intron,fwer)]
	cluster_pval[,padj:=p.adjust(fwer,'fdr')]

	sig=foreach(i=fdr,.combine='rbind')%do%{
		sig=sum(unlist(cluster_pval[,padj<i]),na.rm=TRUE)
		data.table(fdr=i,sig=sig)
	}

	fwrite(sig,sprintf('%s/%s.sig.txt',out_dir,out_prefix),sep='\t')

}

# Main: 
main(in_fn,out_dir,out_prefix)

