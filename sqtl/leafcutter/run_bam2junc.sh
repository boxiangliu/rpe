# convert bam to junction:
bam_dir=$1
out_dir=$2
mkdir -p $out_dir

if [[ -f $out_dir/juncfiles.txt ]]; then rm $out_dir/juncfiles.txt; fi
n=0

for bamfile in `ls $bam_dir/*_dedup.rg.bam`; do
	n=$((n+1))
	if [[ n -ge 10 ]]; then 
		wait
		n=0
	fi
	juncfile=$out_dir/$(basename $bamfile | cut -d'_' -f1,1).junc
	echo Converting $bamfile to $juncfile
	sh /srv/persistent/bliu2/tools/leafcutter/scripts/bam2junc.sh $bamfile $juncfile &
done
wait

ls $out_dir/*.junc > $out_dir/juncfiles.txt
