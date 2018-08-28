# find_motif.sh finds sequences that are recognized as motif 
# matches (based on the motif PWM). 


#!!!! This analysis is not complete yet. 
out_dir=/srv/persistent/bliu2/rpe/processed_data/motif_enrichment/find_motif/
mkdir -p $out_dir

motif_fn=../processed_data/motif_enrichment/motif_enrichment/galactose/homerResults/motif2.motif

perl /srv/persistent/bliu2/tools/HOMER/bin/findMotifs.pl \
$in_dir/glucose.target.fa fasta \
$out_dir/glucose/ \
-fasta $in_dir/glucose.background.fa \
-find $motif_fn 
