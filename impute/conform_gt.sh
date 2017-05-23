java=/srv/gsfs0/software/java/jre1.8.0_66/bin/java
conform_gt=/srv/gsfs0/projects/montgomery/bliu2/tools/conform-gt.24May16.cee.jar
reference_dir=/srv/gsfs0/projects/montgomery/bliu2/shared/beagle/
gt=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/plink2vcf/rpe.merged.missing5e-2.autosome.fill.vcf
out_dir=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/impute

# for i in $(seq 1 22);do
for i in 22;do
# i=22 # test on chr22
$java -Xmx32g -jar $conform_gt chrom=$i ref=$reference_dir/chr$i.1kg.phase3.v5a.vcf.gz strict=false gt=$gt out=$out_dir/rpe.conform_gt.1kgp3v5.chr$i > $out_dir/rpe.conform_gt.1kgp3v5.chr$i.log
done
