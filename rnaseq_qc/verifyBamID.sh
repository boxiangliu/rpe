export out_dir=../processed_data/rnaseq_qc/verifyBamID/
mkdir -p $out_dir

extract_chr1_and_convert_to_GRCh37(){
	in_fn=$1
	out_fn=$2
	samtools view -h $in_fn chr1 | sed 's/chr1/1/' | samtools view -hb - -o $out_fn
	samtools index $out_fn
}
export -f extract_chr1_and_convert_to_GRCh37

rehead_vcf(){
	in_fn=$1
	old_to_new=$2
	out_fn=$3
	bcftools reheader -s $old_to_new -o $out_fn $in_fn 
}

export -f rehead_vcf 

run_verifyBamID(){
	bam=$1
	sample=$(basename $bam)
	extract_chr1_and_convert_to_GRCh37 $bam $out_dir/$sample
	verifyBamID --vcf $out_dir/rpe.imputed.chr1.all_filters.vcf.gz \
	--bam $out_dir/$sample \
	--out $out_dir/${sample}_out \
	--verbose \
	--noPhoneHome \
	--ignoreRG \
	--best
}

export -f run_verifyBamID

rehead_vcf ../data/genotype/filt/rpe.imputed.chr1.all_filters.vcf.gz \
../data/meta/dna2rna.txt \
$out_dir/rpe.imputed.chr1.all_filters.vcf.gz

parallel -j10 run_verifyBamID {} ::: $(ls ../data/rnaseq/bam/*/*.bam)