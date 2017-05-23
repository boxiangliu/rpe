java=/srv/gsfs0/software/java/jre1.8.0_66/bin/java
beagle=/srv/gsfs0/software/beagle/4.1/beagle.12Oct15.b2c.jar
reference_dir=/srv/gsfs0/projects/montgomery/bliu2/shared/beagle/
gt=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/plink2vcf/rpe.merged.missing5e-2.autosome.fill.vcf
out_dir=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/impute

# for i in $(seq 1 22);do
for i in 22;do
# i=22 # test on chr22
$java -Xmx32g -jar $beagle nthreads=10 chrom=$i ref=$reference_dir/chr$i.1kg.phase3.v5a.vcf.gz map=$reference_dir/plink.chr$i.GRCh37.map impute=true gt=$gt out=$out_dir/rpe.beagle.1kgp3v5.chr$i > $out_dir/rpe.beagle.1kgp3v5.chr$i.log
done