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
python rasqual/calc_gcc.py /mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf exon > ../processed_data/rasqual/gcc.exon.txt
bash rasqual/htseq.sh
Rscript rasqual/htseq.merge.R
Rscript rasqual/prepare_input_files.R

# Run rasqual:
parallel -j10 bash rasqual/rasqual.sh ../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} ../processed_data/rasqual/expression/glucose.expression.bin ../processed_data/rasqual/expression/glucose.size_factors_gc.bin ../data/genotype/asvcf/glucose/rpe.imputed.chr22.all_filters.vcf.new.gz chr22 '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.log ::: {1..396}


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
