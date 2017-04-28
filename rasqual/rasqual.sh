# run RASQUAL
# bash rasqual.sh <input parameter file> <line number> Y.bin K.bin VCF X.bin

param_file=$1
line_num=$2
Y=$3
K=$4
vcf_file=$5
mode=$6 # total, ase, or joint
out_dir=$7
X=${8:-null}


param=($(cat $1 | sed "${line_num}q;d"))
gene_id=${param[0]}
gene_name=${param[1]}
region=${param[2]}
n_rsnp=${param[3]}
n_fsnp=${param[4]}
exon_start_positions=${param[5]}
exon_end_positions=${param[6]}
feat_id=$(grep $gene_id -n ../processed_data/rasqual/expression/glucose.expression.txt | cut -d":" -f1,1) # The line number corresponding to the gene_id.
window_size=2000000
n_sample=24
echo id: $gene_id 
echo name: $gene_name 
echo region: $region
echo reference snps: $n_rsnp
echo feature snps: $n_fsnp
echo feature id: $feat_id



if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi


if [[ $X == 'null' ]]; then 
	echo INFO - running without covariates.
	if [[ -e $out_dir/${gene_id}_${gene_name}.txt ]]; then 
		echo "$out_dir/${gene_id}_${gene_name}.txt exist! skipping..."
	else
		if [[ $mode == 'joint' ]]; then 
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			-y $Y \
			-k $K \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		elif [[ $mode == 'total' ]]; then 
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			--population-only \
			-y $Y \
			-k $K \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		elif [[ $mode == 'ase' ]]; then
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			--as-only \
			-y $Y \
			-k $K \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		else 
			echo '--mode must be either total, ase, or joint!'
		fi
	fi

else
	echo INFO - running with covariates $X.
	if [[ -e $out_dir/${gene_id}_${gene_name}.txt ]]; then 
		echo "$out_dir/${gene_id}_${gene_name}.txt exist! skipping..."
	else
		if [[ $mode == 'joint' ]]; then 
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			-y $Y -k $K -x $X \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		elif [[ $mode == 'total' ]]; then 
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			--population-only \
			-y $Y -k $K -x $X \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		elif [[ $mode == 'ase' ]]; then
			tabix $vcf_file $region | \
			/srv/persistent/bliu2/tools/rasqual/bin/rasqual \
			--as-only \
			-y $Y -k $K -x $X \
			-n $n_sample -j $feat_id -l $n_rsnp -m $n_fsnp \
			-s $exon_start_positions -e $exon_end_positions \
			--imputation-quality 0.8 --imputation-quality-fsnp 0.8 \
			--minor-allele-frequency 0.05 --hardy-weinberg-pvalue 0.0 \
			--minor-allele-frequency-fsnp 0.05 \
			--cis-window-size $window_size \
			-f ${gene_id}_${gene_name} --n_threads 1 \
			--force -v --genotype-dosage > $out_dir/${gene_id}_${gene_name}.txt
		else 
			echo '--mode must be either total, ase, or joint!'
		fi
	fi
fi 