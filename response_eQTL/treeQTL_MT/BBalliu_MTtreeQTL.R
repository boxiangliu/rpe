#!/bin/R
# r v 3.3.1 
# Brunilda Balliu
# July 25th 2017
# Script to run sinlge and multi-tissue eQTL analysis for the RPE project - Adapted from scripts I got from Bosh
# In durga

library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)

# Directories:
rasqual_dir_glucose='/srv/persistent/bliu2/rpe/processed_data/rasqual/output/glucose/joint/'
out_dir_glucose='/srv/persistent/bliu2/rpe/processed_data/rasqual/output/glucose/treeQTL/'
rasqual_dir_galactose='/srv/persistent/bliu2/rpe/processed_data/rasqual/output/galactose/joint/'
out_dir_galactose='/srv/persistent/bliu2/rpe/processed_data/rasqual/output/galactose/treeQTL/'
out_dir='/srv/persistent/bballiu/RPE'
setwd(out_dir)

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

MTtreeQTL=function(meqtl_glucose, meqtl_galactose, snp_map, gene_map, level1=0.05,level2=0.05,level3=0.05,eSNP=FALSE, out_dir,cis_dist=1e6) {
  fwrite(meqtl_glucose,sprintf('%s/mEQTL_glucose_out_cis.txt',out_dir),sep='\t')
  # on.exit(unlink(sprintf('%s/mEQTL_glucose_out_cis.txt',out_dir)))
  
  fwrite(meqtl_galactose,sprintf('%s/mEQTL_galactose_out_cis.txt',out_dir),sep='\t')
  # on.exit(unlink(sprintf('%s/mEQTL_galactose_out_cis.txt',out_dir)))
  
  snps_by_tissue=data.frame(snp_name=snp_map$snpid, glucose=rep(1,nrow(snp_map)), galactose=rep(1,nrow(snp_map)))
  genes_by_tissue=data.frame(gene_name=gene_map$geneid, glucose=rep(1,nrow(gene_map)), galactose=rep(1,nrow(gene_map)))
  
  if (eSNP){
    print('INFO - level1 is eSNP')
    
    # Get number of genes nearby to each SNP:
    print('INFO - calculating number of tests per SNP...')
    n_tests_per_SNP=get_n_tests_per_SNP(snp_map,gene_map,nearby=TRUE,dist=1e+06)
    
    # Get list of eSNPs:
    print('INFO - getting eSNPs...')
    eSNPs=get_eSNPs_multi_tissue(genes_by_tissue=genes_by_tissue, snps_by_tissue=snps_by_tissue, n_tests_per_SNP=n_tests_per_SNP, m_eqtl_out_dir=out_dir, tissue_names=c("galactose","glucose"), level1 = 0.05, level2 = 0.05, level3 = 0.05)    
    return(eSNPs)
    
  } else {
    print('INFO - level1 is eGene')
    # Get list of eGenes:
    print('INFO - getting eGenes and eAssociations...')
    eGenes=get_eGenes_multi_tissue(genes_by_tissue = genes_by_tissue, snps_by_tissue = snps_by_tissue, snp_map = snp_map, gene_map = gene_map, nearby = TRUE, dist = 1e+06, m_eqtl_out_dir=out_dir,  tissue_names=c("galactose","glucose"), level1 = 0.05, level2 = 0.05, level3 = 0.05)
    return(eGenes)
  }
}

# Read RASQUAL result: 
tmp_file_glucose=tempfile()
system(sprintf("cat %s/*/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",rasqual_dir_glucose,tmp_file_glucose))
joint_glucose=fread(tmp_file_glucose,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file_glucose)

tmp_file_galactose=tempfile()
system(sprintf("cat %s/*/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",rasqual_dir_galactose,tmp_file_galactose))
joint_galactose=fread(tmp_file_galactose,col.names=c('fid','sid','log10qval','chisq','pi'))
file.remove(tmp_file_galactose)

# Take the ensembl ID as fid:
joint_glucose[,fid:=str_split_fixed(fid,'_',2)[,1]]
joint_galactose[,fid:=str_split_fixed(fid,'_',2)[,1]]

# Calculate p-value from chisq stat: 
joint_glucose[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]
joint_galactose[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]

# Reorder based on p-value: 
setorder(joint_glucose,pval)
setorder(joint_galactose,pval)

# Multiple hypothesis correction using TreeQTL:
gene_map=get_gene_map(unique(joint_glucose$fid))
snp_map=get_snp_map(joint_glucose$sid)

# Create MatrixEQTL format output 
meqtl_glucose=joint_glucose[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
meqtl_galactose=joint_galactose[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]

# Run multi-tissue (i.e. multi-treatment) treeQTL and save results 
# This code only works if the directory "out_dir" only contains the meqtl_glucose and meqtl_galactose results.
eGenes=MTtreeQTL(meqtl_glucose, meqtl_galactose, snp_map, gene_map, level1=0.05,level2=0.05,level3=0.05, eSNP=FALSE,out_dir,cis_dist=1e6)
fwrite(eGenes,sprintf('%s/eGenesMT.txt',out_dir),sep='\t')
table(eGenes[,2:3]) # number of treatments each gene is active on

# Run single tissue treeQTL and save results
temp_glucose=treeQTL(meqtl_glucose,snp_map,gene_map,level1=.05,level2=.05,eSNP=FALSE)
eGenes_glucose=temp_glucose[[1]]
eAssocs_glucose=temp_glucose[[2]]
setorder(eGenes_glucose,fam_p)
setorder(eAssocs_glucose,BBFDR)
fwrite(eGenes_glucose,sprintf('%s/eGenes_glucose.txt',out_dir),sep='\t')
fwrite(eAssocs_glucose,sprintf('%s/eAssocs_glucose.txt',out_dir),sep='\t')

temp_galactose=treeQTL(meqtl_galactose,snp_map,gene_map,level1=.05,level2=.05,eSNP=FALSE)
eGenes_galactose=temp_galactose[[1]]
eAssocs_galactose=temp_galactose[[2]]
setorder(eGenes_galactose,fam_p)
setorder(eAssocs_galactose,BBFDR)
fwrite(eGenes_galactose,sprintf('%s/eGenes_galactose.txt',out_dir),sep='\t')
fwrite(eAssocs_galactose,sprintf('%s/eAssocs_galactose.txt',out_dir),sep='\t')

# Naive validation approach (see report)
GeneNames=unique(meqtl_glucose$gene)
JointDist=matrix(c(table(GeneNames %in% eGenes_glucose$family,GeneNames %in% eGenes_galactose$family)), nrow = 2, dimnames = list("Age 70" = c("No", "Yes"), "Age 80" = c("No", "Yes")))
mcnemar.test(x=JointDist,correct = T)

# Two-step FDR cutoff approach from Yoav Gilad's paper (see report)
# Validate glucose eGenes on the galactose data
meqtl_gluINgal=meqtl_galactose[meqtl_galactose$gene %in% eGenes_glucose$family,]
gene_map_gluINgal=data.frame(gene_map[gene_map$geneid %in% eGenes_glucose$family, ])
temp_gluINgal=treeQTL(meqtl = meqtl_gluINgal, snp_map = snp_map, gene_map = gene_map_gluINgal,  level1=.5,level2=.5,eSNP=FALSE)
eGenes_gluINgal=temp_gluINgal[[1]]
eAssocs_gluINgal=temp_gluINgal[[2]]
setorder(eGenes_gluINgal,fam_p)
setorder(eAssocs_gluINgal,BBFDR)
fwrite(eGenes_gluINgal,sprintf('%s/eGenes_gluINgalFDR50.txt',out_dir),sep='\t')
fwrite(eAssocs_gluINgal,sprintf('%s/eAssocs_gluINgalFDR50.txt',out_dir),sep='\t')

meqtl_galINglu=meqtl_glucose[meqtl_glucose$gene %in% eGenes_galactose$family,]
gene_map_galINglu=data.frame(gene_map[gene_map$geneid %in% eGenes_galactose$family, ])
temp_galINglu=treeQTL(meqtl = meqtl_galINglu, snp_map = snp_map, gene_map = gene_map_galINglu, level1=.5,level2=.5,eSNP=FALSE)
eGenes_galINglu=temp_galINglu[[1]]
eAssocs_galINglu=temp_galINglu[[2]]
setorder(eGenes_galINglu,fam_p)
setorder(eAssocs_galINglu,BBFDR)
fwrite(eGenes_galINglu,sprintf('%s/eGenes_galINgluFDR50.txt',out_dir),sep='\t')
fwrite(eAssocs_galINglu,sprintf('%s/eAssocs_galINgluFDR50.txt',out_dir),sep='\t')

eGenesGluGal=rbind(c(nrow(eGenes_glucose),nrow(eGenes_gluINgal)), c(nrow(eGenes_galactose),nrow(eGenes_galINglu)))
eGenesGluGal[,2]/eGenesGluGal[,1]

eAssocsGluGal=rbind(c(nrow(eAssocs_glucose),nrow(eAssocs_gluINgal)), c(nrow(eAssocs_galactose),nrow(eAssocs_galINglu)))

