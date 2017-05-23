out_dir=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/rnaseq/sorted
processed_data=/srv/gsfs0/projects/montgomery/bliu2/rpe/processed_data/
n=0
for input in `ls $processed_data/rnaseq/bam/set*/*Aligned.out.bam`; do
	n=n+1
	base=$(basename $input)
	echo $base
	samtools sort -o $out_dir/${base/out.bam/out.sorted.bam} $input &
	if [[ $n -ge 10 ]];then 
		wait 
		n=0
	fi 
done 