library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)

# Default:
rasqual_dir='../processed_data/rasqual/output/glucose/joint/'
out_dir='../processed_data/rasqual/output/glucose/treeQTL/'
level1 = 0.05
level2 = 0.05

args=commandArgs(T)
if (length(args)>0){
	rasqual_dir=args[1]
	out_dir=args[2]
	level1=as.numeric(args[3])
	level2=as.numeric(args[4])
}

if (!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)

# Functions:
get_gene_map=function(gene_ids){
	gene_ids=unique(gene_ids)
	gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')
	gencode[,geneid:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
	gene_map=unique(gencode[V3=='gene',.(geneid,chr=V1,left=V4,right=V5)])
	gene_map=gene_map[geneid%in%gene_ids]
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
	on.exit(unlink(sprintf('%s/mEQTL_out_cis.txt',tmp_dir)))

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
		eAssocs=get_eAssociations(eSNPs,n_tests_per_SNP,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=TRUE)

		return(list(eSNPs=eSNPs,eAssocs=eAssocs))

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
		eAssocs=get_eAssociations(eGenes, n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir), sprintf("%s/eAssoc_cis.txt",tmp_dir),by_snp=FALSE)
		
		return(list(eGenes=eGenes,eAssocs=eAssocs))
	}
}


# Read RASQUAL result: 
tmp_file=tempfile()
system(sprintf("cat %s/*/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",rasqual_dir,tmp_file))
joint=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file)


# Take the ensembl ID as fid:
joint[,fid:=str_split_fixed(fid,'_',2)[,1]]


# Calculate p-value from chisq stat: 
joint[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


# Reorder based on p-value: 
setorder(joint,pval)


# Multiple hypothesis correction using TreeQTL:
gene_map=get_gene_map(unique(joint$fid))
snp_map=get_snp_map(joint$sid)

meqtl=joint[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
temp=treeQTL(meqtl,snp_map,gene_map,level1=level1,level2=level2,eSNP=FALSE)
eGenes=temp[[1]]
eAssocs=temp[[2]]
setorder(eGenes,fam_p)
setorder(eAssocs,BBFDR)

# Save results:
fwrite(eGenes,sprintf('%s/eGenes.txt',out_dir),sep='\t')
fwrite(eAssocs,sprintf('%s/eAssocs.txt',out_dir),sep='\t')