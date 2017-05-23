in_dir=../processed_data/rasqual/input/
expr_dir=../processed_data/rasqual/expression/
geno_dir=../data/genotype/asvcf/galactose_nodup/
out_dir=../processed_data/rasqual/output/galactose/joint/
log_dir=../logs/rasqual/
cov_dir=../processed_data/select_covariates/select_covariates/covariates/
start=$1
end=$2
i=$3


for j in `seq $start $end`; do
	if [[ ! -d $out_dir/chr$i/ ]]; then mkdir -p $out_dir/chr$i/; fi
	if [[ ! -d $log_dir/chr$i/ ]]; then mkdir -p $log_dir/chr$i/; fi
	bash rasqual/scg4/rasqual.sh \
		$in_dir/rasqual.input.chr$i.filt.txt \
		$j \
		$expr_dir/galactose.expression.bin \
		$expr_dir/galactose.size_factors_gc.bin \
		$geno_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
		joint \
		$out_dir/chr$i/ \
		$cov_dir/gal_cov.9.bin
done
