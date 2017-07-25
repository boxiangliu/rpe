in_dir='../data/genotype/asvcf/glucose_nodup/'
out_dir='../processed_data/matrixEQTL/genotype/'
[[ ! -d $out_dir ]] && mkdir -p $out_dir

for i in `seq 22`; do
bcftools view -q 0.05 -Ou $in_dir/rpe.imputed.chr$i.all_filters.vcf.new.gz | \
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' \
-o $out_dir/chr$i.txt
done



cat $out_dir/chr{1..22}.txt | sed -re 's/\[[0-9]+\]//g' -e 's/:DS//g' | \
awk 'BEGIN{FS="\t";OFS="\t"}{out=$1"_"$2"_"$3"_"$4"_b37";for (i=5;i<=NF;i++){out=out"\t"$i}; print out}' | \
sed 's/# CHROM_POS_REF_ALT_b37/snp_id/' | \
awk 'BEGIN{FS="\t";OFS="\t";count=0}{if ($1=="snp_id") {count=count+1; if (count==1) {print $0;}} else {print $0;}}' > $out_dir/genotype.txt

cat $out_dir/chr{1..22}.txt | sed '/^#/d' | \
awk 'BEGIN{FS="\t";OFS="\t";print "snp_id","chr","pos"}{snp_id=$1"_"$2"_"$3"_"$4"_b37";print snp_id,$1,$2}' > $out_dir/snpsloc.txt