# Make ASVCF: 
echo "counting ASE..."
parallel -j6 bash /srv/persistent/bliu2/tools/rasqual/src/ASVCF/createASVCF.sh rasqual/glucose.bam_list.txt ../data/genotype/rpe.imputed.chr{}.dr2.vcf.gz ::: {1..22}

