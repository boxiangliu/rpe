out_dir=../processed_data/genotype_qc/detect_duplication/
mkdir -p $out_dir

bcftools query -H \
-f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' \
../data/genotype/filt/rpe.imputed.chr1.all_filters.vcf.gz > \
$out_dir/chr1_dosage.txt

Rscript genotype_qc/detect_duplication.R