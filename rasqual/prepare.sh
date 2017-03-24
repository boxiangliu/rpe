# Copy BAM files from scg4:
while read line; do
	scp bliu2@carmack.stanford.edu:$line /srv/persistent/bliu2/rpe/data/rnaseq/
	scp bliu2@carmack.stanford.edu:${line/.bam/.bai/} /srv/persistent/bliu2/rpe/data/rnaseq/
done < rasqual/bam_list.scg4.txt


# Move file to appropriate directories: 
cd /srv/persistent/bliu2/rpe/data/rnaseq
mv *lucose*ba* glucose/
mv *alactose*ba* galactose/

# Link genotype files: 
cp /users/xli6/data/xin/rpe/imputation/beagle/rpe.imputed.chr*.vcf.gz /srv/persistent/bliu2/rpe/data/genotype/orig

# Filter VCF using GTEx standard pipeline:
tail -n +2 ../data/meta/dna2rna.txt > /srv/scratch/bliu2/rpe.genotype.reheader.txt
tail -n +2 ../data/meta/dna2rna.txt | sort -k2 | cut -f1 > /srv/scratch/bliu2/rpe.genotype.sample_order.txt

for i in {1..22}; do
	tabix -p vcf ../data/genotype/orig/rpe.imputed.chr$i.vcf.gz
	bcftools view -m2 -M2 -Ov -S /srv/scratch/bliu2/rpe.genotype.sample_order.txt ../data/genotype/orig/rpe.imputed.chr$i.vcf.gz | bcftools reheader -s /srv/scratch/bliu2/rpe.genotype.reheader.txt | awk 'BEGIN{OFS="\t"}{if ($1 !~ /#/) {print "chr"$0 } else {print $0}}' | awk 'BEGIN{OFS="\t"}{if ($1 !~ /#/) {if (length($5)<=51 && length($4)<=51) {print $0}} else {print $0}}' | bcftools filter -e "AR2<0.8" -Oz -o /srv/persistent/bliu2/rpe/data/genotype/filt/rpe.imputed.chr$i.all_filters.vcf.gz
	tabix -p vcf ../data/genotype/filt/rpe.imputed.chr$i.all_filters.vcf.gz
done
