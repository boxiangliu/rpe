in_dir=../processed_data/rasqual/input/
expr_dir=../processed_data/rasqual/expression/
geno_dir=../data/genotype/asvcf/glucose_nodup/
out_dir=../processed_data/rasqual/output/glucose/joint/
log_dir=../logs/rasqual/
cov_dir=../processed_data/select_covariates/select_covariates/covariates/

# for i in {1..22}; do 
for i in {15..22}; do 
	echo INFO - chr$i

	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
	echo INFO - $n_genes genes.

	if [[ ! -d $out_dir/chr$i/ ]]; then mkdir -p $out_dir/chr$i/; fi
	if [[ ! -d $log_dir/chr$i/ ]]; then mkdir -p $log_dir/chr$i/; fi

	# parallel -j15 bash rasqual/rasqual.sh \
	# 	$in_dir/rasqual.input.chr$i.filt.txt \
	# 	{} \
	# 	$expr_dir/glucose.expression.bin \
	# 	$expr_dir/glucose.size_factors_gc.bin \
	# 	$geno_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
	# 	joint \
	# 	$out_dir/chr$i/ \
	# 	$cov_dir/glu_cov.8.bin '2>' $log_dir/chr$i/rasqual.chr22.{}.log ::: `seq $n_genes`

	n=0
	for j in `seq $n_genes`; do
		n=$((n+1))
		if [[ n -gt 15 ]]; then wait; n=0; fi
		bash rasqual/rasqual.sh \
			$in_dir/rasqual.input.chr$i.filt.txt \
			$j \
			$expr_dir/glucose.expression.bin \
			$expr_dir/glucose.size_factors_gc.bin \
			$geno_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
			joint \
			$out_dir/chr$i/ \
			$cov_dir/glu_cov.8.bin 2> $log_dir/chr$i/rasqual.chr22.$j.log &
	done
done