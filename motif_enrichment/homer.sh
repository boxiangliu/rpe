in_dir=/srv/persistent/bliu2/rpe/processed_data/motif_enrichment/get_sequence/
out_dir=/srv/persistent/bliu2/rpe/processed_data/motif_enrichment/motif_enrichment/
mkdir -p $out_dir

for condition in glucose galactose; do
	perl /srv/persistent/bliu2/tools/HOMER/bin/findMotifs.pl \
	$in_dir/$condition.target.fa fasta \
	$out_dir/$condition/ \
	-fasta $in_dir/$condition.background.fa
done