in_dir='../data/genotype/asvcf/glucose_nodup/'
out_dir='../processed_data/matrixEQTL/genotype/'
[[ ! -d $out_dir ]] && mkdir -p $out_dir

for i in `seq 22`; do
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' \
$in_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz \
-o $out_dir/chr$i.txt
done 