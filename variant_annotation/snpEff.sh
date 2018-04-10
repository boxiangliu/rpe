mkdir -p ../data/genotype/filt/snpEff/
out_dir=../processed_data/variant_annotation/snpEff/
mkdir -p $out_dir

# Functions: 
extract_anno(){

	in_fn=$1
	anno=$2
	out_fn=$3

	echo INFO - input: $in_fn
	echo INFO - anno: $anno
	echo INFO - output: $out_fn

	java -Xmx4g \
	-jar /srv/persistent/bliu2/tools/snpEff/SnpSift.jar \
	filter "ANN[ANY].EFFECT has '$anno'" \
	$in_fn | \
	grep -v ^# | \
	awk -v anno=$anno 'BEGIN{OFS="\t"}{print $1,$2,$2,anno}' | \
	bgzip > $out_fn

}

export -f extract_anno


# Annotate using SnpEff:
parallel -j22 \
java -Xmx4g \
-jar /srv/persistent/bliu2/tools/snpEff/snpEff.jar \
GRCh37.75 \
../data/genotype/filt/rpe.imputed.chr{}.all_filters.vcf.gz '|' \
bgzip '>' \
../data/genotype/filt/snpEff/chr{}.snpEff.vcf.gz ::: {1..22}

# Annote with ClinVar:
parallel -j22 \
java -Xmx4g \
-jar /srv/persistent/bliu2/tools/snpEff/SnpSift.jar \
annotate ../data/clinvar/clinvar.vcf.gz \
../data/genotype/filt/snpEff/chr{}.snpEff.vcf.gz '|' \
bgzip '>' \
../data/genotype/filt/snpEff/chr{}.snpEff.clinvar.vcf.gz ::: {1..22}

# Remove intermediate files: 
rm ../data/genotype/filt/snpEff/chr*.snpEff.vcf.gz

# Concatenate:
zcat ../data/genotype/filt/snpEff/chr{1..22}.snpEff.clinvar.vcf.gz | \
bgzip > ../data/genotype/filt/snpEff/all.snpEff.clinvar.vcf.gz

# Extract all annotations using SnpSift:
parallel -j10 \
extract_anno ../data/genotype/filt/snpEff/all.snpEff.clinvar.vcf.gz \
{} \
$out_dir/{}.bed.gz \
::: downstream_gene_variant exon_variant intron_variant missense_variant \
splice_acceptor_variant splice_donor_variant splice_region_variant \
synonymous_variant upstream_gene_variant 3_prime_UTR_variant 5_prime_UTR_variant

# Concatenate:
zcat $out_dir/*_variant.bed.gz | bgzip > $out_dir/all.bed.gz
