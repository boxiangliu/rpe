log_dir=../logs/
size=100
for i in `seq 22`; do
	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
	for start in `seq 1 $size $n_genes`; do
		end=$((start+size-1))
		if [[ $end -gt $n_genes ]]; then end=$n_genes; fi
		echo INFO - chr$i:$start-$end 
		qsub -N ras.chr$i.$start-$end -l h_vmem=8G -cwd -o $log_dir/rasqual/qsub/rasqual.chr$i.$start-$end.log -e $log_dir/rasqual/qsub/rasqual.chr$i.$start-$end.err rasqual/scg4/rasqual.driver.sh $start $end $i
	done
done 