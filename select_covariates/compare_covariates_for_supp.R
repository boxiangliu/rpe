library(data.table)
library(stringr)
library(cowplot)
library(dplyr)
library(TreeQTL)

in_dir_glu='../processed_data/select_covariates/select_covariates/glucose/'
in_dir_gal='../processed_data/select_covariates/select_covariates/galactose/'
fig_dir = '../figures/select_covariates/compare_covariates_for_supp/'

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

treeQTL_eSNP=function(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,tmp_dir='./',cis_dist=1e6){
	fwrite(meqtl,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sep='\t')

	# Get number of genes nearby to each SNP:
	print('INFO - calculating number of tests per SNP...')
	n_tests_per_SNP=get_n_tests_per_SNP(snp_map,gene_map,nearby=TRUE,dist=cis_dist)

	# Get list of eSNPs:
	print('INFO - getting eSNPs...')
	eSNPs=get_eSNPs(n_tests_per_SNP, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)

	unlink(sprintf('%s/mEQTL_out_cis.txt',tmp_dir))
	return(eSNPs)
}

treeQTL_eGenes=function(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,tmp_dir='./',cis_dist=1e6){
	fwrite(meqtl,sprintf('%s/mEQTL_out_cis.txt',tmp_dir),sep='\t')

	# Get number of SNPs nearby to each gene:
	print('INFO - calculating number of tests per gene...')
	n_tests_per_gene=get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = cis_dist)


	# Get list of eGenes:
	print('INFO - getting eGenes...')
	eGenes=get_eGenes(n_tests_per_gene, sprintf('%s/mEQTL_out_cis.txt',tmp_dir),level1=level1,level2=level2)


	unlink(sprintf('%s/mEQTL_out_cis.txt',tmp_dir))
	return(eGenes)
}


#------------- Glucose -------------
# Read RASQUAL result: 
tmp_file=tempfile()
joint_ls=list()
for (n in 1:8){
	print(sprintf('INFO - %s covariates..',n))
	# Read RASQUAL result:
	system(sprintf("cat %s/chr22/cov%s/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",in_dir_glu,n,tmp_file))
	joint=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))


	# Take the ensembl ID as fid:
	joint[,fid:=str_split_fixed(fid,'_',2)[,1]]


	# Calculate p-value from chisq stat: 
	joint[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


	# Reorder based on p-value: 
	setorder(joint,pval)
	joint_ls[[n]]=joint
}
file.remove(tmp_file)


# Multiple hypothesis correction using TreeQTL:
gene_map=get_gene_map()
joint_adj_ls=list()

for (i in 1:length(joint_ls)){
	joint=joint_ls[[i]]
	meqtl=joint[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
	snp_map=get_snp_map(meqtl$SNP)
	joint_adj=treeQTL(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,eSNP=TRUE)
	setorder(joint_adj,BBFDR)
	joint_adj_ls[[i]]=joint_adj
}


# Count the number of significant associations:
cov=c('sex',paste0('genoPC',1:3),paste0('sva',1:4))
sig=data.frame()
for (i in 1:length(joint_ls)){
	sig=rbind(sig,data.frame(cov=cov[i],n=sum(joint_adj_ls[[i]]$BBFDR<0.05)))
}
p1=ggplot(sig,aes(cov,n))+
	geom_point(size=3)+
	geom_line(aes(group=1))+
	ggtitle('Glucose')+
	theme(axis.text.x=element_text(angle=45,vjust=0.5))+
	xlab('') +
	ylab('# Significant eAssociations\n(FDR < 0.05)')+
	scale_x_discrete(breaks = c('sex',paste0('genoPC',1:3),paste0('sva',1:4)),labels = c('sex',paste0('geno PC',1:3),paste0('SV',1:4)))

#------------- Galactose ------------------
# Read RASQUAL result: 
tmp_file=tempfile()
joint_ls=list()
for (n in 1:9){
	print(sprintf('INFO - %s covariates..',n))
	# Read RASQUAL result:
	system(sprintf("cat %s/chr22/cov%s/*.txt | awk '{print $1,$2,$10,$11,$12}' > %s",in_dir_gal,n,tmp_file))
	joint=fread(tmp_file,col.names=c('fid','sid','log10qval','chisq','pi'))


	# Take the ensembl ID as fid:
	joint[,fid:=str_split_fixed(fid,'_',2)[,1]]


	# Calculate p-value from chisq stat: 
	joint[,pval:=pchisq(q=chisq,df=1,lower.tail=F)]


	# Reorder based on p-value: 
	setorder(joint,pval)
	joint_ls[[n]]=joint
}
file.remove(tmp_file)


# Multiple hypothesis correction using TreeQTL:
gene_map=get_gene_map()
joint_adj_ls=list()

for (i in 1:length(joint_ls)){
	joint=joint_ls[[i]]
	meqtl=joint[,list(SNP=sid,gene=fid,beta=-1,`t-stat`=-1,`p-value`=pval,FDR=-1)]
	snp_map=get_snp_map(meqtl$SNP)
	joint_adj=treeQTL(meqtl,snp_map,gene_map,level1=0.05,level2=0.05,eSNP=TRUE)
	setorder(joint_adj,BBFDR)
	joint_adj_ls[[i]]=joint_adj
}


# Count the number of significant associations:
cov=c('sex',paste0('genoPC',1:3),paste0('sva',1:5))
sig=data.frame()
for (i in 1:length(joint_ls)){
	sig=rbind(sig,data.frame(cov=cov[i],n=sum(joint_adj_ls[[i]]$BBFDR<0.05)))
}
p2=ggplot(sig,aes(cov,n))+
	geom_point(size=3)+
	geom_line(aes(group=1))+
	ggtitle('Galactose')+
	theme(axis.text.x=element_text(angle=45,vjust=0.5))+
	xlab('') +
	ylab('# Significant eAssociations\n(FDR < 0.05)')+
	scale_x_discrete(breaks = c('sex',paste0('genoPC',1:3),paste0('sva',1:5)),labels = c('sex',paste0('geno PC',1:3),paste0('SV',1:5)))


p = plot_grid(p1,p2,labels = c('A','B'), align = 'h', nrow = 1)
fig_fn = sprintf('%s/eAssociation_vs_cov_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8, base_height = 4)


p = plot_grid(p1,p2,labels = c('a','b'), align = 'h', nrow = 1)
fig_fn = sprintf('%s/eAssociation_vs_cov_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 8, base_height = 4)
