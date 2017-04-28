in_dir=../data/genotype/asvcf/glucose_sid
for f in $(ls $in_dir/rpe.imputed.chr*.all_filters.vcf.new.gz); do
	echo $f
	bcftools view -s ^021011 $f -Oz > ${f/_sid/_nodup}
	tabix -p vcf ${f/_sid/_nodup}
done 


in_dir=../data/genotype/asvcf/galactose_sid
for f in $(ls $in_dir/rpe.imputed.chr*.all_filters.vcf.new.gz); do
	echo $f
	bcftools view -s ^021011 $f -Oz > ${f/_sid/_nodup}
	tabix -p vcf ${f/_sid/_nodup}
done 