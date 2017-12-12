# Calculate splice levels:


# link Glucose and Galactose junction files into the same directory:
mkdir -p /srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/both/

for fn in 020206-2 020310 020311 021010 021512 021611 030411 032411 041212-2944 041212-9319 051211 052711-1 052711-2 070810 071709 072910 080410 081209 081309 081508 081909 082609 091109; do
	ln ../data/rnaseq/leafcutter/glucose/$fn.junc ../data/rnaseq/leafcutter/both/${fn}_glucose.junc
	ln ../data/rnaseq/leafcutter/galactose/$fn.junc ../data/rnaseq/leafcutter/both/${fn}_galactose.junc
done 

# Count splice levels:
ls /srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/both/*.junc > /srv/persistent/bliu2/rpe/data/rnaseq/leafcutter/both/juncfiles.txt
mkdir -p ../data/rnaseq/leafcutter/both/cluster/
python /srv/persistent/bliu2/tools/leafcutter/clustering/leafcutter_cluster.py \
	-j ../data/rnaseq/leafcutter/both/juncfiles.txt \
	-r ../data/rnaseq/leafcutter/both/cluster/ \
	-o diff_splice