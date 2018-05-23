out_dir='../processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/'
hapmap_fn='/srv/persistent/bliu2/shared/ldscore/w_hm3.snplist'
amd_gwas_fn='../data/gwas/Fritsche_2015_AdvancedAMD.v2.txt'
myopia_gwas_fn='../data/gwas/23andme_myopia.prepared.txt.gz'

mkdir -p $out_dir

Rscript disease_enrichment/ld_score_regression/order_amd_gwas.R

~/tools/ldsc/munge_sumstats.py \
	--sumstats $amd_gwas_fn \
	--out $out_dir/amd \
	--merge-alleles $hapmap_fn \
	--N-cas-col Ncases \
	--N-con-col Ncontrols \
	--p GC.Pvalue \
	--a1 effect_allele \
	--a2 non_effect_allele \
	--snp Marker \
	--a1-inc

~/tools/ldsc/munge_sumstats.py \
	--sumstats $myopia_gwas_fn \
	--out $out_dir/myopia \
	--merge-alleles $hapmap_fn \
	--p pvalue \
	--a1 ref \
	--a2 alt \
	--snp rsid \
	--signed-sumstats beta,0 \
	--N-cas 106086 \
	--N-con 85757
