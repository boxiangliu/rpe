out_dir=../processed_data/sqtl/optimal_covariate_sva/test_covariate/
mkdir -p $out_dir

# Remove some introns that causes FastQTL to error:
mkdir -p ../data/rnaseq/leafcutter/galactose/cluster_clean
zcat ../data/rnaseq/leafcutter/galactose/cluster/sqtl_perind.counts.gz.qqnorm_chr1.gz | \
grep -v "chr1:44069204:44070574:clu_15512" | grep -v "chr1:26140661:26142039:clu_15247" | \
bgzip > ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr1.gz
tabix -p bed ../data/rnaseq/leafcutter/galactose/cluster_clean/sqtl_perind.counts.gz.qqnorm_chr1.gz

mkdir -p ../data/rnaseq/leafcutter/glucose/cluster_clean
cp ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr1.gz ../data/rnaseq/leafcutter/glucose/cluster_clean
cp ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr1.gz.tbi ../data/rnaseq/leafcutter/glucose/cluster_clean

# Functions:
run_FastQTL(){

n=$1 # number of covariates
i=$2 # chromosome 
cov=$3 # covariate file
bed_dir=$4 # bed directory
vcf_dir=$5
out_dir=$6
echo $out_dir
mkdir -p $out_dir

if [[ n -gt 0 ]]; then

head -n $((n+1)) $cov > $out_dir/temp_covariate.chr$i.$n.txt
/users/zappala/software/fastqtl/bin/fastQTL \
--seed 42 \
--maf-threshold 0.05 \
--vcf $vcf_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
--bed $bed_dir/sqtl_perind.counts.gz.qqnorm_chr$i.gz \
--cov $out_dir/temp_covariate.chr$i.$n.txt \
--region chr$i \
--window 1e5 \
--permute 1000 \
--out $out_dir/chr$i.permute.$n.txt.gz \
--log $out_dir/chr$i.permute.$n.log

else

/users/zappala/software/fastqtl/bin/fastQTL \
--seed 42 \
--maf-threshold 0.05 \
--vcf $vcf_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
--bed $bed_dir/sqtl_perind.counts.gz.qqnorm_chr$i.gz \
--region chr$i \
--window 1e5 \
--permute 1000 \
--out $out_dir/chr$i.permute.$n.txt.gz \
--log $out_dir/chr$i.permute.$n.log

fi
}

# Main
export -f run_FastQTL
parallel -j15 run_FastQTL \
{1} \
1 \
../processed_data/optimal_covariate_sva/combine_covariates/{2}_covariate.txt \
../data/rnaseq/leafcutter/{2}/cluster_clean/ \
../data/genotype/asvcf/{2}_nodup/ \
$out_dir/{2}/ \
::: {0..6} ::: glucose galactose

# Remove temporary directories:
rm -r ../data/rnaseq/leafcutter/{glucose,galactose}/cluster_clean/
