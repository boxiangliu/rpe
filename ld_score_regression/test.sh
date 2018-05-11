python /users/bliu2/tools/ldsc/munge_sumstats.py \
--sumstats /srv/persistent/bliu2/rpe/data/ldsc/BMI/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt \
--merge-alleles /srv/persistent/bliu2/rpe/data/ldsc/w_hm3.snplist \
--out /srv/persistent/bliu2/rpe/data/ldsc/BMI/BMI \
--a1-inc

python /users/bliu2/tools/ldsc/ldsc.py \
  --h2 /srv/persistent/bliu2/rpe/data/ldsc/BMI/BMI.sumstats.gz \
  --ref-ld-chr /srv/persistent/bliu2/rpe/data/ldsc/1000G_Phase3_baselineLD_ldscores/baselineLD. \
  --frqfile-chr /srv/persistent/bliu2/rpe/data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
  --w-ld-chr /srv/persistent/bliu2/rpe/data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --overlap-annot \
  --print-coefficients \
  --print-delete-vals \
  --out /srv/persistent/bliu2/rpe/data/ldsc/BMI/BMI.baselineLD

Rscript /users/bliu2/tools/ldsc/ContinuousAnnotations/quantile_h2g.r \
/srv/persistent/bliu2/rpe/data/ldsc/1000G_Phase3_baselineLD_ldscores/baselineLD.MAF_Adj_Predicted_Allele_Age.q5.M \
/srv/persistent/bliu2/rpe/data/ldsc/BMI/BMI.baselineLD \
/srv/persistent/bliu2/rpe/data/ldsc/BMI/BMI.baselineLD.Allele_Age.q5.txt
