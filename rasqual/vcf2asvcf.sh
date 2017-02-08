# Make ASVCF: 
echo "counting ASE..."
# parallel -j6 bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh rasqual/glucose.bam_list.txt ../data/genotype/filt/rpe.imputed.chr{}.all_filters.vcf.gz ::: {1..22}
parallel -j6 bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh rasqual/galactose.bam_list.txt ../data/genotype/filt/rpe.imputed.chr{}.all_filters.vcf.gz ::: {1..22}

