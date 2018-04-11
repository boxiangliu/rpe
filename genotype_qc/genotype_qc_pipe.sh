# Post-imputation filtering:
bash genotype_qc/post_imputation_filtering.sh

# Calculate basic stats with bcftools:
bash genotype_qc/bcftools_stats.sh 

# Plot Ts/Tv ratio:
Rscript genotype_qc/plot_tstv.R

# Detect duplicate:
bash genotype_qc/detect_duplication.sh 