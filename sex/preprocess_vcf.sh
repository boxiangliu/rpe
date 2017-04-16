out_dir=../processed_data/sex/preprocess_vcf/
[[ ! -d $out_dir ]] && mkdir -p $out_dir

# Subset to chromosome 1 and X: 
bcftools view -r 1,23 ../data/genotype/before_imputation/rpe.merged.missing5e-2.vcf.gz | sed -r -e 's/^23\t/X\t/' | bgzip > $out_dir/rpe.chr1_and_X.vcf.gz
tabix -p vcf $out_dir/rpe.chr1_and_X.vcf.gz

# Extract genotype dosage:
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $out_dir/rpe.chr1_and_X.vcf.gz | sed -r -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/1/2/g' -e 's/\.\/\./NA/g' -e 's/\[[0-9]{1,2}\]//g' -e 's/:GT//g' -e 's/# //' > $out_dir/rpe.chr1_and_X.tsv