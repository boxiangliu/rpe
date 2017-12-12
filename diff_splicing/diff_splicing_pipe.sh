# Perform differential splicing analysis
# Boxiang Liu
# 2017-12-11

# Calculate splice levels:
bash diff_splicing/preprocess.sh

# Make sample table:
Rscript diff_splicing/make_sample_file.R 


# Do differential splicing:
Rscript diff_splicing/diff_splicing.R
