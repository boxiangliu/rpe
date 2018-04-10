filter(){
i=$1
bcftools view -m2 -M2 --exclude "INFO/AR2<0.8" ../data/genotype/orig/rpe.imputed.chr$i.vcf.gz -Oz -o ../data/genotype/filt/rpe.imputed.chr$i.all_filters.vcf.gz
tabix -p vcf ../data/genotype/filt/rpe.imputed.chr$i.all_filters.vcf.gz
}
export -f filter
parallel -j22 filter {} ::: {1..22}