mkdir -p ../processed_data/sqtl/fastQTL/{nominal,permutation}/{glucose,galactose}

###########
# Glucose #
###########
# Nominal:
parallel -j15 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--region chr{} \
--window 1e5 \
--out ../processed_data/sqtl/fastQTL/nominal/glucose/chr{}.nominal.txt.gz ::: {1..22}

zcat ../processed_data/sqtl/fastQTL/nominal/glucose/chr{1..22}.nominal.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz

# Permutation: 
# Note: try different seeds for failed jobs (seed=7 works for chr20, 102 for chr15, 203 for chr13, 406 for chr5)
parallel -j15 /users/zappala/software/fastqtl/bin/fastQTL \
--seed 42 \
--vcf ../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--region chr{} \
--maf-threshold 0.05 \
--window 1e5 \
--permute 1000 \
--log ../processed_data/sqtl/fastQTL/permutation/glucose/chr{}.permutation.log \
--out ../processed_data/sqtl/fastQTL/permutation/glucose/chr{}.permutation.txt.gz ::: {1..22}

zcat ../processed_data/sqtl/fastQTL/permutation/glucose/chr{1..22}.permutation.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/permutation/glucose/all.permutation.txt.gz

#############
# Galactose #
#############

# Nominal:
parallel -j15 /users/zappala/software/fastqtl/bin/fastQTL \
--vcf ../data/genotype/asvcf/galactose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--region chr{} \
--window 1e5 \
--out ../processed_data/sqtl/fastQTL/nominal/galactose/chr{}.nominal.txt.gz ::: {1..22}

zcat ../processed_data/sqtl/fastQTL/nominal/galactose/chr{1..22}.nominal.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/nominal/galactose/all.nominal.txt.gz



# Permutation: 
# Note: try different seeds for failed jobs (seed=7 works for chr17, 102 for chr16)
parallel -j15 /users/zappala/software/fastqtl/bin/fastQTL \
--seed 42 \
--vcf ../data/genotype/asvcf/galactose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--region chr{} \
--maf-threshold 0.05 \
--window 1e5 \
--permute 1000 \
--log ../processed_data/sqtl/fastQTL/permutation/galactose/chr{}.permutation.log \
--out ../processed_data/sqtl/fastQTL/permutation/galactose/chr{}.permutation.txt.gz ::: 1 15

mkdir -p ../data/rnaseq/leafcutter/galactose/cluster_clean/
zcat ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr1.gz | \
grep -v "chr1:26135280:26135586:clu_15246" |\
bgzip > ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr1.gz
tabix -p bed ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr1.gz

zcat ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr15.gz |\
grep -v "chr15:83656136:83657770:clu_6022" |\
grep -v "chr15:89450607:89456478:clu_6063" |\
bgzip > ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr15.gz 
tabix -f -p bed ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr15.gz

parallel -j15 /users/zappala/software/fastqtl/bin/fastQTL \
--seed 42 \
--vcf ../data/genotype/asvcf/galactose_nodup/rpe.imputed.chr{}.all_filters.vcf.new.gz \
--bed ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr{}.gz \
--region chr{} \
--maf-threshold 0.05 \
--window 1e5 \
--permute 1000 \
--log ../processed_data/sqtl/fastQTL/permutation/galactose/chr{}.permutation.log \
--out ../processed_data/sqtl/fastQTL/permutation/galactose/chr{}.permutation.txt.gz ::: 15


zcat ../processed_data/sqtl/fastQTL/permutation/galactose/chr{1..22}.permutation.txt.gz | \
gzip > ../processed_data/sqtl/fastQTL/permutation/galactose/all.permutation.txt.gz

