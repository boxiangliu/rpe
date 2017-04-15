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


#-------- RASQUAL --------------
# Setup: 
mkdir rasqual ../processed_data/rasqual


# Prepare data for RASQUAL:
bash rasqual/addreadgroup.sh
bash rasqual/prepare.sh
bash rasqual/vcf2asvcf.sh
bash rasqual/change_sid.sh
python rasqual/calc_gcc.py /mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf exon > ../processed_data/rasqual/gcc.exon.txt
bash rasqual/htseq.sh
Rscript rasqual/htseq.merge.R
Rscript rasqual/prepare_input_files.R

# Test RASQUAL (only 50 genes in joint, total, and ase mode):
bash rasqual/test/test_rasqual.sh
bash rasqual/test/test_rasqual_chr22.sh

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


#-------- genotype PC -------------
# Setup: 
mkdir -p genotype_pc ../processed_data/genotype_pc ../figures/genotype_pc


# Pick random 1000 Genome samples:
Rscript genotype_pc/1kg_samples.R


# Extract dosage from VCF files: 
bash genotype_pc/preprocess_vcf.sh


# Do PCA on genotype dosage: 
Rscript genotype_pc/genotype_pc.R


# 