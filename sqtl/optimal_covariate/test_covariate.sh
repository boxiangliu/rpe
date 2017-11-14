geno_pc=$1
splice_pc=$2
peer=$3
vcf=$4
bed=$5
out_dir=$6
mkdir -p $out_dir 

# peer=../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.PCs
# vcf=../data/genotype/asvcf/glucose_nodup/rpe.imputed.chr22.all_filters.vcf.new.gz
# bed=../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind.counts.gz.qqnorm_chr22.gz

echo INFO - genotype PC = $geno_pc
echo INFO - splicing PC = $splice_pc


echo INFO - making covariates

Rscript sqtl/utils/combine_covariates.R \
	--genotype_pc=../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv \
	--peer=$peer \
	--gender=../processed_data/sex/sex/gender.tsv \
	--output=$out_dir/covariates-${geno_pc}_geno_pc-${splice_pc}_expr_pc.tsv \
	--gender_coding=letter \
	--num_geno_pc=$geno_pc \
	--num_peer_factor=$splice_pc \
	--row_and_colnames=TRUE


echo INFO - running FastQTL
/users/zappala/software/fastqtl/bin/fastQTL \
--vcf $vcf \
--bed $bed \
--cov $out_dir/covariates-${geno_pc}_geno_pc-${splice_pc}_expr_pc.tsv \
--region chr22 \
--window 1e5 \
--out $out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.txt.gz ::: {1..22}



echo INFO - counting significant loci
Rscript sqtl/optimal_covariate/calc_FDR.R \
$out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.txt.gz \
$out_dir/chr22.nominal.geno_pc-$geno_pc.splice_pc-$splice_pc.sig.txt \
$geno_pc $splice_pc
