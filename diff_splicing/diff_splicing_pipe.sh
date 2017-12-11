# Perform differential splicing analysis
# Boxiang Liu
# 2017-12-11

# Make sample table:
Rscript diff_splicing/make_sample_file.R 

# Do differential splicing:
../scripts/leafcutter_ds.R --num_threads 4 --exons_file=../leafcutter/data/gencode19_exons.txt.gz ../example_data/testYRIvsEU_perind_numers.counts.gz ../example_data/test_diff_intron.txt

