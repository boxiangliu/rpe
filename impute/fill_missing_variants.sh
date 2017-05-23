awk 'BEGIN{OFS="\t"}
{if (NF==1 || $1=="#CHROM") {
	print $0}
else if ($5==".") {
	$5=$4;
	print $0;
}}' /srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/plink2vcf/rpe.merged.missing5e-2.autosome.vcf > /srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/plink2vcf/rpe.merged.missing5e-2.autosome.fill.vcf
