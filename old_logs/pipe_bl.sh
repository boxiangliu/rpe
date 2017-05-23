#### initialization:
scripts=/srv/gsfs0/projects/montgomery/bliu2/rpe/scripts
processed_data=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data

#### 160819: 
#### convert PLINK to VCF format: 
# setup: 
mkdir $scripts/plink2vcf
mkdir $processed_data/plink2vcf

# back up 17-26/plink2:
cp -r /srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/17-26/plink2 /srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/17-26/plink2_bak


# manually change the name of 17-26/plink2/PLINK_040816_0234/plink2.ped
# 2       052711 -> 02     052711.2
# 2       052711 -> 03     052711.1


# merge three arrays:
echo "/srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/9-16/plink/PLINK_040816_0240/plink.ped /srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/9-16/plink/PLINK_040816_0240/plink.map" > $processed_data/plink2vcf/merge_list.txt
echo "/srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/17-26/plink2/PLINK_040816_0234/plink2.ped /srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/17-26/plink2/PLINK_040816_0234/plink2.map" >> $processed_data/plink2vcf/merge_list.txt
plink --file /srv/gsfs0/projects/montgomery/bliu2/rpe/raw/hfRPE_Genotyping/plink/1-8/plink/PLINK_040816_0225/plink \
--merge-list $processed_data/plink2vcf/merge_list.txt \
--recode --out $processed_data/plink2vcf/rpe.merged


# filter SNPs with missingness rate > 5%: 
plink --file $processed_data/plink2vcf/rpe.merged --geno 0.05 --recode --out $processed_data/plink2vcf/rpe.merged.missing5e-2


# convert plink file to vcf: 
plink --file $processed_data/plink2vcf/rpe.merged.missing5e-2 --recode vcf-iid --out $processed_data/plink2vcf/rpe.merged.missing5e-2

#### impute
# subset to autosomes: 
bash impute/subset_autosome.sh 

# fill in missing alternative genotype:
bash impute/fill_missing_variants.sh

# imputation with BEAGLE:
bash $scripts/impute/download_beagle_ref.sh

#### calculate RPKMs
# set up
mkdir $processed_data/rnaseq
mkdir $processed_data/rnaseq/sorted
mkdir $scripts/rnaseq/
ln -s /srv/gsfs0/projects/montgomery/nsabell/rpe/bam $processed_data/rnaseq



# sort bam files: 
bash $scripts/rnaseq/sort.sh
# index bam files 
# make sample sheeet
# run rnaseqc 


#### calculate peer factors 
#### calculate genotype PCs
#### call eQTL
#### 