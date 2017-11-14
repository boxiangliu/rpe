# Create directories:
mkdir -p ../processed_data/sqtl/ ../figures/sqtl/


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


# Select optimal combination of covariates:
parallel -j10 --joblog ../processed_data/sqtl/optimal_covariate/test_covariate/log \
bash sqtl/optimal_covariate/test_covariate.sh \
{1} {2} ../processed_data/sqtl/optimal_covariate/test_covariate/ \
::: 2 3 4 ::: 1 2 3 4 5 6 7 8 9 10

cat ../processed_data/sqtl/optimal_covariate/test_covariate/chr22.nominal.*.sig.txt \
> ../processed_data/sqtl/optimal_covariate/test_covariate/sig.txt

Rscript sqtl/optimal_covariate/compare_covariate.R


# Add ancestry PC and sex:
Rscript sqtl/utils/combine_covariates.R \
	--genotype_pc=../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv \
	--peer=../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.PCs \
	--gender=../processed_data/sex/sex/gender.tsv \
	--output=../processed_data/sqtl/covariate/covariates-2_geno_pc-4_splice_pc.tsv \
	--gender_coding=letter \
	--num_geno_pc=2 \
	--num_peer_factor=4 \
	--row_and_colnames=TRUE


# sQTL mapping with fastQTL in nominal mode:
parallel -j10 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--cov ../processed_data/sqtl/covariate/covariates-2_geno_pc-4_splice_pc.tsv \
--region chr{} \
--window 1e5 \
--out ../processed_data/sqtl/fastQTL/nominal/chr{}.nominal.txt.gz ::: {1..22}


zcat ../processed_data/sqtl/fastQTL/nominal/chr{1..22}.nominal.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz



# sQTL mapping with fastQTL in permutation mode:
zcat ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr16.gz | \
sed '/chr16:29461597:29463430:clu_4930/d' | \
bgzip > ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr16.2.gz # remove a line that triggers GSL domain error

mv ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr16.2.gz \
../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr16.gz

tabix -p bed ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr16.gz

parallel -j1 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--cov ../processed_data/sqtl/covariate/covariates-2_geno_pc-4_splice_pc.tsv \
--region chr{} \
--window 1e5 \
--permute 1000 10000 \
--out ../processed_data/sqtl/fastQTL/permutation/chr{}.permutation.txt.gz ::: 16

zcat ../processed_data/sqtl/fastQTL/permutation/chr{1..22}.permutation.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz


# Plot p-values: 
Rscript sqtl/fastQTL/plot_pvalue.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../figures/sqtl/fastQTL/plot_pvalue/


# Plot sQTL vs distance to TSS:
Rscript sqtl/fastQTL/plot_sqtl_vs_distance.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../figures/sqtl/fastQTL/plot_sqtl_vs_distance/


# Count sQTLs: 
Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/ \
fastqtl.nominal

Rscript sqtl/fastQTL/count_sqtl.R \
../processed_data/sqtl/fastQTL/permutation/all.permutation.txt.gz \
../processed_data/sqtl/fastQTL/count_sqtl/ \
fastqtl.permutation