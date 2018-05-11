# Combine RPKMs:
ls ../data/rnaseq/rnaseqc/*/genes.rpkm.gct > gct_paths.txt

Rscript rnaseq/combine_GCTs.R \
gct_paths.txt \
../data/rnaseq/rpkm/rpe.rnaseqc_rpkm

rm gct_paths.txt

# Combine read counts:
ls ../data/rnaseq/rnaseqc/*/*/*.metrics.tmp.txt.intronReport.txt > count_paths.txt

Rscript rnaseq/combine_readCounts.R \
count_paths.txt \
../data/rnaseq/count/rpe_rnaseqc_matrix.txt

rm count_paths.txt