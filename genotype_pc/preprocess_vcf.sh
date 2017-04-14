mkdir -p ../processed_data/genotype_pc/preprocess_vcf/rpe/
mkdir -p ../processed_data/genotype_pc/preprocess_vcf/1kg/

# RPE VCFs: 
parallel -j11 'bcftools query -H -f "%ID[\t%DS]\n" ../data/genotype/asvcf/glucose_sid/rpe.imputed.chr{}.all_filters.vcf.new.gz > ../processed_data/genotype_pc/preprocess_vcf/chr{}.tsv' ::: {1..22}

# 1KG VCFs:
parallel -j11 "bcftools view -S ../processed_data/genotype_pc/1kg_samples/4_sample_each_pop.txt -m2 -M2 -Ou /srv/persistent/bliu2/shared/1000genomes/phase3v5a/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT\_b37' -Ou | bcftools query -H -f '%ID[\t%GT]\n' | sed -r -e 's/0\\|0/0/g' -e 's/0\\|1/1/g' -e 's/1\\|0/1/g' -e 's/1\\|1/2/g' > ../processed_data/genotype_pc/preprocess_vcf/1kg/chr{}.tsv" ::: {1..22}

