# Create directories:
mkdir -p ../processed_data/sqtl/ ../figures/sqtl/

#--------- Preprocess ------------#
#---------- Glucose --------------#
# Convert BAM to junction files:
bash sqtl/leafcutter/run_bam2junc.sh \
/srv/persistent/bliu2/rpe/data/rnaseq/bam/glucose/ \
/srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/glucose/


# Intron clustering:
sed -i '/021011/d' ../data/rnaseq/leafcutter/glucose/juncfiles.txt # remove duplicate
mkdir ../data/rnaseq/leafcutter/glucose/cluster/
python /srv/persistent/bliu2/tools/leafcutter/clustering/leafcutter_cluster.py \
	-j ../data/rnaseq/leafcutter/glucose/juncfiles.txt \
	-r ../data/rnaseq/leafcutter/glucose/cluster/ \
	-o sqtl


# Prepare FastQTL input:
python /srv/persistent/bliu2/tools/leafcutter/scripts/prepare_phenotype_table.py \
	../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz \
	-p 10


# Change GRCh37 to hg19 coordinates:
parallel -j1 \
cat ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{} "|" \
python sqtl/utils/b37_to_hg19.py ">" ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 ::: {1..22}

parallel -j1 \
mv ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 \
../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{} ::: {1..22}


# Bgzip and index the bed files:
bash ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz_prepare.sh # bgzip and index


#---------------- Galactose -----------------#
# Convert BAM to junction files:
bash sqtl/leafcutter/run_bam2junc.sh \
/srv/persistent/bliu2/rpe/data/rnaseq/bam/galactose/ \
/srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/galactose/


# Intron clustering:
sed -i '/021011/d' ../data/rnaseq/leafcutter/galactose/juncfiles.txt # remove duplicate
mkdir ../data/rnaseq/leafcutter/galactose/cluster/
python /srv/persistent/bliu2/tools/leafcutter/clustering/leafcutter_cluster.py \
	-j ../data/rnaseq/leafcutter/galactose/juncfiles.txt \
	-r ../data/rnaseq/leafcutter/galactose/cluster/ \
	-o sqtl


# Prepare FastQTL input:
python /srv/persistent/bliu2/tools/leafcutter/scripts/prepare_phenotype_table.py \
	../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz \
	-p 10


# Change GRCh37 to hg19 coordinates:
parallel -j1 \
cat ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{} "|" \
python sqtl/utils/b37_to_hg19.py ">" ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 ::: {1..22}

parallel -j1 \
mv ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}_2 \
../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{} ::: {1..22}


# Bgzip and index the bed files:
bash ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz_prepare.sh # bgzip and index


#--------- Mapping sQTLs for both conditions -----------#
# Select covariates:
Rscript sqtl/optimal_covariate_sva/combine_covariates.R
bash sqtl/optimal_covariate_sva/test_covariate.sh
Rscript sqtl/optimal_covariate_sva/calc_FDR.R


# sQTL mapping with fastQTL:
bash sqtl/fastQTL/run_fastQTL.sh


# Quality control:
Rscript sqtl/fastQTL/quality_control.R

# Adjust p-value and assign intron to genes:
Rscript sqtl/fastQTL/adjust_pvalue.R

# Pick interesting examples:
Rscript sqtl/fastQTL/sQTL_examples.R


# Count sQTLs: 
Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/glucose/ \
fastqtl.nominal

Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/glucose/ \
fastqtl.permutation

Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/nominal/galactose/all.nominal.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/galactose/ \
fastqtl.nominal


Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/galactose/ \
fastqtl.permutation
