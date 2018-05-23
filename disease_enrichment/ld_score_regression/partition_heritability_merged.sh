unset DISPLAY XAUTHORITY

export amd_sumstats_fn='../processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/amd.sumstats.gz'
export myopia_sumstats_fn='../processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/myopia.sumstats.gz'
export out_dir='../processed_data/disease_enrichment/ld_score_regression/partition_heritability_merged/'
export baseline_annotation_dir='/srv/persistent/bliu2/shared/ldscore/baseline/'
export tissue_specific_annotation_dir='../processed_data/disease_enrichment/ld_score_regression/ldscore_merged/'
export weight_dir='/srv/persistent/bliu2/shared/ldscore/weights_hm3_no_hla/'
export frq_dir='/srv/persistent/bliu2/shared/ldscore/1000G_frq/'

mkdir -p $out_dir


partition_heritability(){
mkdir -p $2
echo INFO - in_dir: $1
echo INFO - out_dir: $2
echo INFO - trait: $3
echo INFO - gwas: $4
if [[ -e "$2/$3.merged.results" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi 

python ~/tools/ldsc/ldsc.py \
--h2 $4 \
--ref-ld-chr $1/merged.,$baseline_annotation_dir/baseline. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $2/$3.merged \
--print-coefficients
}

export -f partition_heritability

partition_heritability_nobaseline(){

mkdir -p $2
echo INFO - in_dir: $1
echo INFO - out_dir: $2
echo INFO - trait: $3
echo INFO - gwas: $4
if [[ -e "$2/$3.merged.nobaseline.results" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi 

python ~/tools/ldsc/ldsc.py \
--h2 $4 \
--ref-ld-chr $1/merged. \
--w-ld-chr $weight_dir/weights. \
--overlap-annot \
--frqfile-chr $frq_dir/1000G.mac5eur. \
--out $2/$3.merged.nobaseline \
--print-coefficients
}

export -f partition_heritability_nobaseline

parallel -j15 partition_heritability $tissue_specific_annotation_dir/{} $out_dir/{} amd $amd_sumstats_fn ::: 4sd top200 top500 top1000 top2000 top4000
parallel -j15 partition_heritability_nobaseline $tissue_specific_annotation_dir/{} $out_dir/{} amd $amd_sumstats_fn ::: 4sd top200 top500 top1000 top2000 top4000
parallel -j15 partition_heritability $tissue_specific_annotation_dir/{} $out_dir/{} myopia $myopia_sumstats_fn ::: 4sd top200 top500 top1000 top2000 top4000
parallel -j15 partition_heritability_nobaseline $tissue_specific_annotation_dir/{} $out_dir/{} myopia $myopia_sumstats_fn ::: 4sd top200 top500 top1000 top2000 top4000
