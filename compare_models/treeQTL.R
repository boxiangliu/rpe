library(data.table)
library(TreeQTL)
library(stringr)

in_fn='../processed_data/compare_models/gtex_eqtl/all_tissues.chr22.lowest_pval.txt'
out_dir='../processed_data/compare_models/gtex_eqtl/'
if (!dir.exists(out_dir)) dir.create(out_dir)

# Read Cis-eQTL and coerce to Matrix eQTL format:
cis=fread(in_fn,col.names=c('gene','SNP','dist','p-value','beta','se','tissue'))
cis[,c("t-stat","FDR"):=-1]
cis=cis[,.(SNP,gene,beta,`t-stat`,`p-value`,FDR)]


# Order by p-value:
setorder(cis,`p-value`)


# Out to tmp file: 
fwrite(cis,'mEQTL_out_cis.txt',sep='\t')


# Distance used to define nearby region:
dist=1e6


# Get SNP map: 
snp_map=data.table(snpid=unique(cis$SNP))
snp_map[,chr:=str_split_fixed(snpid,'_',5)[,1]]
snp_map[,chr:=paste0('chr',chr)]
snp_map[,pos:=as.integer(str_split_fixed(snpid,'_',5)[,2])]


# Get gene map: 
gencode=fread('/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf')
gencode[,geneid:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
gene_map=unique(gencode[V3=='gene',.(geneid,chr=V1,left=V4,right=V5)])


# Get number of genes nearby to each SNP:
n_tests_per_SNP=get_n_tests_per_SNP(snp_map, gene_map, nearby = TRUE, dist = dist)


# Get list of eSNPs:
eSNPs=get_eSNPs(n_tests_per_SNP, "mEQTL_out_cis.txt")


# Generate txt output file with full set of eAssociations:
get_eAssociations(eSNPs, n_tests_per_SNP, "mEQTL_out_cis.txt", sprintf("%s/eAssoc_cis_by_snp.txt",out_dir),by_snp=TRUE)


# Get number of SNPs nearby to each gene:
n_tests_per_gene=get_n_tests_per_gene(snp_map, gene_map, nearby = TRUE, dist = dist)


# Get list of eGenes:
eGenes=get_eGenes(n_tests_per_gene, "mEQTL_out_cis.txt")


# Generate txt output file with full set of eAssociations:
get_eAssociations(eGenes, n_tests_per_gene, "mEQTL_out_cis.txt", sprintf("%s/eAssoc_cis_by_gene.txt",out_dir),by_snp=FALSE)


# Remove intermediate files: 
unlink('mEQTL_out_cis.txt')