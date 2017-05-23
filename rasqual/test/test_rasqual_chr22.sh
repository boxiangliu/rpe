#------------- all chr22, using VCF files with chr_pos_ref_alt_b37 as RSID, no covariates ---------------
mkdir -p ../processed_data/rasqual/test_chr22/
# Run rasqual (total read count mode):
echo -e '[total mode]\nstart time' > ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose_sid/rpe.imputed.chr22.all_filters.vcf.new.gz \
total ../processed_data/rasqual/test_chr22/total/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.total.log ::: {1..396}
echo 'end time' >> ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt


# Run rasqual (ase mode):
echo -e '[ase mode]\nstart time' >> ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose_sid/rpe.imputed.chr22.all_filters.vcf.new.gz \
ase ../processed_data/rasqual/test_chr22/ase/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.log ::: {1..396}
echo 'end time' >> ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt


# Run rasqual (joint mode):
echo -e '[joint mode]\nstart time' >> ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose_sid/rpe.imputed.chr22.all_filters.vcf.new.gz \
joint ../processed_data/rasqual/test_chr22/joint/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.log ::: {1..396}
echo 'end time' >> ../processed_data/rasqual/test_chr22/profile.txt
date >> ../processed_data/rasqual/test_chr22/profile.txt
