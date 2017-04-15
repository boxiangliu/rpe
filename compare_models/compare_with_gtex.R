library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)

fig_dir='../figures/compare_models/compare_with_gtex/'
if (!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE)

# Functions:
get_gene_map=function(){
	gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')
	gencode[,geneid:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
	gene_map=unique(gencode[V3=='gene',.(geneid,chr=V1,left=V4,right=V5)])
	return(gene_map)
}

get_snp_map=function(snps){
	snp_map=data.table(snpid=unique(snps))
	snp_map[,chr:=str_split_fixed(snpid,'_',5)[,1]]
	snp_map[,chr:=paste0('chr',chr)]
	snp_map[,pos:=as.integer(str_split_fixed(snpid,'_',5)[,2])]
	return(snp_map)
}


treeQTL=function(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,eSNP=TRUE,tmp_dir='./',cis_dist=1e6){
	fwrite(meqtl,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sep='\t')

	if (eSNP){
		print('INFO - level1 is eSNP')

		# Get number of genes nearby to each SNP:
		print('INFO - calculating number of tests per SNP...')
		n_tests_per_SNP=get_n_tests_per_SNP(snp_map,gene_map,nearby=TRUE,dist=cis_dist)

		# Get list of eSNPs:
		print('INFO - getting eSNPs...')
		eSNPs=get_eSNPs(n_tests_per_SNP, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)


		# Generate txt output file with full set of eAssociations:
		print('INFO - getting eAssociations...')
		get_eAssociations(eSNPs,n_tests_per_SNP,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=TRUE)

	} else {
		print('INFO - level1 is eGene')

		# Get number of SNPs nearby to each gene:
		print('INFO - calculating number of tests per gene...')
		n_tests_per_gene=get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = cis_dist)


		# Get list of eGenes:
		print('INFO - getting eGenes...')
		eGenes=get_eGenes(n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)


		# Generate txt output file with full set of eAssociations:
		print('INFO - getting eAssociations...')
		get_eAssociations(eGenes, n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir), sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=FALSE)
	}
	res=fread(sprintf("%s/eAssoc_cis.txt",tmp_dir))
	print('INFO - removing intermediate files...')
	unlink(sprintf('%s/mEQTL_out_cis.txt',tmp_dir))
	unlink(sprintf('%s/eAssoc_cis.txt',tmp_dir))
	return(res)
}


# Read RASQUAL result: 
tmp_file=tempfile()
system(sprintf("cat ../processed_data/rasqual/test_chr22/total/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
bi=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
system(sprintf("cat ../processed_data/rasqual/test_chr22/ase/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
as=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
system(sprintf("cat ../processed_data/rasqual/test_chr22/joint/chr22/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",tmp_file))
joint=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file)


# Take the ensembl ID as fid:
bi[,fid:=str_split_fixed(fid,'_',2)[,1]]
as[,fid:=str_split_fixed(fid,'_',2)[,1]]
joint[,fid:=str_split_fixed(fid,'_',2)[,1]]


# Calculate p-value from chisq stat: 
bi[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
as[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
joint[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


# Reorder based on p-value: 
setorder(bi,pval)
setorder(as,pval)
setorder(joint,pval)


# Multiple hypothesis correction using TreeQTL:
gene_map=get_gene_map()

meqtl=bi[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
snp_map=get_snp_map(meqtl$SNP)
bi_adj=treeQTL(meqtl,snp_map,gene_map,level1=35,level2=35,eSNP=TRUE)
setorder(bi_adj,BBFDR)

meqtl=as[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
snp_map=get_snp_map(meqtl$SNP)
as_adj=treeQTL(meqtl,snp_map,gene_map,level1=35,level2=35,eSNP=TRUE)
setorder(as_adj,BBFDR)

meqtl=joint[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
snp_map=get_snp_map(meqtl$SNP)
joint_adj=treeQTL(meqtl,snp_map,gene_map,level1=35,level2=35,eSNP=TRUE)
setorder(joint_adj,BBFDR)


# Load GTEx eQTL (most significant p-value across tissue and corrected with TreeQTL):
gtex=fread('../processed_data/compare_models/gtex_eqtl/eAssoc_cis_by_snp.txt')
true=gtex[BBFDR<0.05,]


# Add ID column: 
gtex[,id:=paste(gene,SNP,sep='_')]
bi_adj[,id:=paste(gene,SNP,sep='_')]
as_adj[,id:=paste(gene,SNP,sep='_')]
joint_adj[,id:=paste(gene,SNP,sep='_')]



# Generate sensitivity vs FDR data: 
FDR=seq(0,0.1,1e-4)
sens=matrix(NA,nrow=length(FDR),ncol=3)
colnames(sens)=c('bi','as','joint')


for (i in 1:length(FDR)){
	sub=bi_adj[BBFDR<FDR[i],id]
	sens[i,'bi']=length(intersect(gtex$id,sub))# takes a while

	sub=as_adj[BBFDR<FDR[i],id]
	sens[i,'as']=length(intersect(gtex$id,sub))# takes a while

	sub=joint_adj[BBFDR<FDR[i],id]
	sens[i,'joint']=length(intersect(gtex$id,sub))# takes a while
}
sens=sens/length(gtex$id)


# Plot:
pdf(sprintf('%s/sensitivity_vs_fdr.pdf',fig_dir))
plot(FDR,sens[,'joint'],type='l',col='black',xlab='FDR',ylab='sensitivity')
lines(FDR,sens[,'as'],col='red')
lines(FDR,sens[,'bi'],col='blue')
legend('topleft',col=c('black','red','blue'),lty=1,legend=c('joint','as','bi'))
dev.off()