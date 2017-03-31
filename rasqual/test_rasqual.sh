# Run rasqual (only 50 genes, total read count mode):
echo -e '[total mode]\nstart time' > ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose/rpe.imputed.chr22.all_filters.vcf.new.gz \
chr22 total ../processed_data/rasqual/test/total/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.total.log ::: {1..50}
echo 'end time' >> ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt


# Run rasqual (only 50 genes, ase mode):
echo -e '[ase mode]\nstart time' >> ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose/rpe.imputed.chr22.all_filters.vcf.new.gz \
chr22 ase ../processed_data/rasqual/test/ase/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.log ::: {1..50}
echo 'end time' >> ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt


# Run rasqual (only 50 genes, joint mode):
echo -e '[joint mode]\nstart time' >> ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt
parallel -j10 bash rasqual/rasqual.sh \
../processed_data/rasqual/input/rasqual.input.chr22.filt.txt {} \
../processed_data/rasqual/expression/glucose.expression.bin \
../processed_data/rasqual/expression/glucose.size_factors_gc.bin \
../data/genotype/asvcf/glucose/rpe.imputed.chr22.all_filters.vcf.new.gz \
chr22 joint ../processed_data/rasqual/test/joint/chr22/ '2>' ../logs/rasqual/chr22/rasqual.chr22.{}.log ::: {1..50}
echo 'end time' >> ../processed_data/rasqual/test/profile.txt
date >> ../processed_data/rasqual/test/profile.txt

