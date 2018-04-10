mkdir -p ../data/genotype/orig/qc
parallel -j 22 bcftools stats ../data/genotype/orig/rpe.imputed.chr{}.vcf.gz '>' ../data/genotype/orig/qc/rpe.imputed.chr{}.vcf.gz.stats ::: {1..22}

mkdir -p ../data/genotype/filt/qc
parallel -j 22 bcftools stats ../data/genotype/filt/rpe.imputed.chr{}.all_filters.vcf.gz '>' ../data/genotype/filt/qc/rpe.imputed.chr{}.all_filters.vcf.gz.stats ::: {1..22}
