# RPE project
# Author: Boxiang Liu
# Email: jollier.liu@gmail.com

#--------- RNAseq -----------
# Setup: 
mkdir rnaseq

# Run RNAseQC:
bash rnaseq/rnaseqc.sh /srv/persistent/bliu2/rpe/data/rnaseq/bam/glucose /srv/persistent/bliu2/rpe/data/rnaseq/rnaseqc/glucose/
bash rnaseq/rnaseqc.sh /srv/persistent/bliu2/rpe/data/rnaseq/bam/galactose /srv/persistent/bliu2/rpe/data/rnaseq/rnaseqc/galactose/


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
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/glucose/joint/ ../processed_data/rasqual/output/glucose/treeQTL/ 0.05 0.05 > treeQTL.glucose.log 
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/galactose/joint/ ../processed_data/rasqual/output/galactose/treeQTL/ 0.05 0.05 > treeQTL.galactose.log

Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/glucose/joint/ ../processed_data/rasqual/output/glucose/treeQTL/fdr0.5 0.5 0.5 > treeQTL.glucose.level1-0.5.level2-0.5.log 
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/galactose/joint/ ../processed_data/rasqual/output/galactose/treeQTL/fdr0.5 0.5 0.5 > treeQTL.galactose.level1-0.5.level2-0.5.log

Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/glucose/joint/ ../processed_data/rasqual/output/glucose/treeQTL/fdr1.0 1.0 1.0 > treeQTL.glucose.level1-1.0.level2-1.0.log 
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/galactose/joint/ ../processed_data/rasqual/output/galactose/treeQTL/fdr1.0 1.0 1.0 > treeQTL.galactose.level1-1.0.level2-1.0.log

Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/glucose/joint/ ../processed_data/rasqual/output/glucose/treeQTL/fdr20.0 20.0 20.0 > treeQTL.glucose.level1-20.0.level2-20.0.log 
Rscript rasqual/treeQTL.R ../processed_data/rasqual/output/galactose/joint/ ../processed_data/rasqual/output/galactose/treeQTL/fdr20.0 20.0 20.0 > treeQTL.galactose.level1-20.0.level2-20.0.log
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


#------------ response eQTL -----------
# Response eQTL based on FDR threshold:
Rscript response_eQTL/fdr_threshold.R 

# Response eQTL based on expression difference:
Rscript response_eQTL/residual.R
Rscript response_eQTL/difference.R
Rscript response_eQTL/matrixeQTL_difference.R


#----------- Metasoft ------------
## Run Metasoft (full sample, select eQTL with p<1e-3): 
# Copy HCASMC eQTL data to appropriate location: 
gzip $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt
ln -s $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz $processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY/HCASMC_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz

# Concatenate eQTL result for all tissue, append tissue name, and sort according to gene ID and SNP: 
bash $scripts/160805/concatenate_eqtl_tissues.sh 

# Make input file for Metasoft:
mkdir $processed_data/160805/metasoft_input/
cat $processed_data/160805/v6p_fastQTL_allpairs_FOR_QC_ONLY2/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py > $processed_data/160805/metasoft_input/metasoft_input.txt

# Split metasoft input by chromosome (to reduce memory footprint for next step): 
bash $scripts/160805/split_metasoft_input_by_chr.sh

# Run METASOFT:
mkdir $processed_data/160805/metasoft_output/
parallel -j12 bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input/metasoft_input.{}.txt $processed_data/160805/metasoft_output/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output/metasoft_output.{}.mcmc.log ::: {1..22} X

# merge metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt


## Run Metasoft (subsampled to 52, select eQTL with p<1e-3): 
# copy HCASMC eQTL data to appropriate location:
mkdir $processed_data/160816/subsampling/HCASMC
ln -s $processed_data/160805/hcasmc.eqtl.pc4.peer8.b37.txt.gz $processed_data/160816/subsampling/HCASMC/HCASMC_52.allpairs.txt.gz

# Concatenate eQTL result for all tissue, append tissue name, and sort according to gene ID and SNP: 
bash $scripts/160805/concatenate_eqtl_tissues.subsample.sh 

# Generate Metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52/ 
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py > $processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt

# Split metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt' ::: {1..22} X

# Run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

## Run Metasoft (subsample to 52, select eQTL with p<1e-2)
# generate metasoft input file:
mkdir $processed_data/160805/metasoft_input_subsample_52_p1e-2/
cat /srv/persistent/bliu2/HCASMC_eQTL/processed_data/160816/subsampling/All_Tissues.allpairs.sorted.txt | python $scripts/160805/gen_metasoft_input.py 1e-2 > $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt

# split Metasoft input by chromosome: 
parallel 'grep "_{}_" ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.txt > ../processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt' ::: {1..22} X

# run METASOFT:
mkdir $processed_data/160805/metasoft_output_subsample_52_p1e-2
parallel bash $scripts/160805/metasoft.core.sh $processed_data/160805/metasoft_input_subsample_52_p1e-2/metasoft_input.{}.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.txt $processed_data/160805/metasoft_output_subsample_52_p1e-2/metasoft_output.{}.mcmc.log ::: {1..22} X

# Merge Metasoft output: 
head -n1 $processed_data/160805/metasoft_output/metasoft_output.1.mcmc.txt > $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt
cat $processed_data/160805/metasoft_output/metasoft_output.{1..22}.mcmc.txt | grep -v RSID >> $processed_data/160805/metasoft_output/metasoft_output.1_22.mcmc.txt

# plot heatmap of m-values: 
Rscript $scripts/160805/plot_mvalue_heatmap.R 


#------------ Finemap ---------------
# Setup:
mdkir finemap

# Make Gviz tracks
bash finemap/atacseq_overlap.gviz.R