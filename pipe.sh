# RPE project
# Author: Boxiang Liu
# Email: jollier.liu@gmail.com


#--------- GWAS ATACseq overlap -------------
# Setup: 
mkdir gwas_atacseq_overlap ../figures/gwas_atacseq_overlap


# Download Encode RPE DHS data: 
bash gwas_atacseq_overlap/download.sh


# Overlap GWAS and RPE ATACseq and compare to Encode and Roadmap samples:
Rscript gwas_atacseq_overlap/overlap.gwas_thresholding.R 

#-------- genotype PC -------------
# Setup: 
mkdir -p genotype_pc ../processed_data/genotype_pc ../figures/genotype_pc


# Pick random 1000 Genome samples:
Rscript genotype_pc/1kg_samples.R


# Extract dosage from VCF files: 
bash genotype_pc/preprocess_vcf.sh


# Do PCA on genotype dosage: 
Rscript genotype_pc/genotype_pc.R


#----------- sex ---------------
# Setup: 
mkdir -p sex ../processed_data/sex ../figures/sex


# Extract chromsome X: 
bash sex/preprocess_vcf.sh 


# Determine sex: 
Rscript sex/sex.R


#----------- select covariates -----------
# Setup: 
mkdir -p select_covariates ../processed_data/select_covariates ../figures/select_covariates


# Merge covariates: 
Rscript select_covariates/merge_covariates.R


# Test different combinations of covariates: 
Rscript select_covariates/select_covariates.R


# Compare covariates: 
Rscript select_covariates/compare_covariates.R


#-------- RASQUAL --------------
# Setup: 
mkdir rasqual ../processed_data/rasqual


# Prepare data for RASQUAL:
bash rasqual/addreadgroup.sh
bash rasqual/prepare.sh
bash rasqual/vcf2asvcf.sh
bash rasqual/change_sid.sh
bash rasqual/remove_dup_from_vcf.sh
python rasqual/calc_gcc.py /mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf exon > ../processed_data/rasqual/gcc.exon.txt
bash rasqual/htseq.sh
Rscript rasqual/htseq.merge.R
Rscript rasqual/prepare_input_files.R

# Test RASQUAL (only 50 genes in joint, total, and ase mode):
bash rasqual/test/test_rasqual.sh
bash rasqual/test/test_rasqual_chr22.sh


# Refer to "compare models" and "select covariates"
# Running RASQUAL:
bash rasqual/rasqual.wrapper.sh
bash rasqual/scg4/rasqual.wrapper.sh # on scg4


# Estimate the number of eQTLs with TreeQTL:
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/glucose/joint/ ../processed_data/rasqual/output/glucose/treeQTL/ > treeQTL.glucose.log 
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/galactose/joint/ ../processed_data/rasqual/output/galactose/treeQTL/ > treeQTL.galactose.log


#------- FastQTL ---------------
# Setup: 
mkdir -p fastqtl ../processed_data/fastqtl/{expression,eqtl} ../figures/fastqtl


# Prepare data for FastQTL:
Rscript fastqtl/count2rpkm.R
Rscript fastqtl/normalize.R
Rscript fastqtl/mat2bed.R
bgzip ../processed_data/fastqtl/expression/glucose.bed && tabix -p bed ../processed_data/fastqtl/expression/glucose.bed.gz 
bgzip ../processed_data/fastqtl/expression/galactose.bed && tabix -p bed ../processed_data/fastqtl/expression/galactose.bed.gz


# Run FastQTL:
bash fastqtl/run_fastqtl.sh
Rscript fastqtl/plot_pvalue_distribution.R


# Compare FastQTL and RASQUAL:
Rscript compare_fastQTL_and_RASQUAL.R


#-------- compare models -----------
# Setup: 
mkdir -p compare_models ../processed_data/compare_models ../figures/compare_models


# Define GTEX eQTLs:
bash compare_models/define_gtex_eqtl.sh


# Compare RASQUAL and FastQTL with GTEx known eQTLs: 
Rscript compare_models/compare_with_gtex.R


#------------- MDS --------------
# Setup: 
mkdir -p mds

# Make MDS plots: 
nohup Rscript mds/preprocess.R > preprocess.log & 
Rscript mds/mds.R
