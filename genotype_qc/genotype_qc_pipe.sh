# Post-imputation filtering:
bash genotype_qc/post_imputation_filtering.sh

# Calculate basic stats with bcftools:
bash genotype_qc/bcftools_stats.sh 

# Count number of variants:
Rscript genotype_qc/count_variant.R


# Plot Ts/Tv ratio:
Rscript genotype_qc/plot_tstv.R

# Detect duplicate:
bash genotype_qc/detect_duplication.sh 