# # First round:
# log_dir=../logs/
# size=100
# for i in `seq 22`; do
# 	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
# 	for start in `seq 1 $size $n_genes`; do
# 		end=$((start+size-1))
# 		if [[ $end -gt $n_genes ]]; then end=$n_genes; fi
# 		echo INFO - chr$i:$start-$end 
# 		qsub -N ras.chr$i.$start-$end -l h_vmem=8G -cwd -o $log_dir/rasqual/qsub/rasqual.chr$i.$start-$end.log -e $log_dir/rasqual/qsub/rasqual.chr$i.$start-$end.err rasqual/scg4/rasqual.driver.sh $start $end $i
# 	done
# done 


# Second round:
# log_dir=../logs/
# size=25
# for i in `seq 22`; do
# 	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
# 	for start in `seq 1 $size $n_genes`; do
# 		end=$((start+size-1))
# 		if [[ $end -gt $n_genes ]]; then end=$n_genes; fi
# 		echo INFO - chr$i:$start-$end 
# 		qsub -N ras.chr$i.$start-$end -l h_vmem=8G -cwd -o $log_dir/rasqual/qsub2/rasqual.chr$i.$start-$end.log -e $log_dir/rasqual/qsub2/rasqual.chr$i.$start-$end.err rasqual/scg4/rasqual.driver.sh $start $end $i
# 	done
# done 


# Third round:
# log_dir=../logs/rasqual/qsub3/
# size=5
# for i in `seq 22`; do
# 	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
# 	for start in `seq 1 $size $n_genes`; do
# 		end=$((start+size-1))
# 		if [[ $end -gt $n_genes ]]; then end=$n_genes; fi
# 		echo INFO - chr$i:$start-$end 
# 		qsub -N ras.chr$i.$start-$end -l h_vmem=8G -cwd -o $log_dir/rasqual.chr$i.$start-$end.log -e $log_dir/rasqual.chr$i.$start-$end.err rasqual/scg4/rasqual.driver.sh $start $end $i
# 	done
# done 


# Fourth round:
# log_dir=../logs/rasqual/qsub4/
# size=1
# for i in `seq 22`; do
# 	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
# 	for line_num in `seq 1 $size $n_genes`; do
# 		param=($(cat ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | sed "${line_num}q;d"))
# 		gene_id=${param[0]}
# 		gene_name=${param[1]}
# 		out_dir=../processed_data/rasqual/output/galactose/joint/chr$i/
# 		if [[ -s $out_dir/${gene_id}_${gene_name}.txt ]]; then
# 			echo "$out_dir/${gene_id}_${gene_name}.txt exist! skipping..."
# 		else
# 			echo INFO - chr$i:$line_num 
# 			qsub -N ras.chr$i.$line_num -l h_vmem=8G -cwd -o $log_dir/rasqual.chr$i.$line_num.log -e $log_dir/rasqual.chr$i.$line_num.err rasqual/scg4/rasqual.driver.sh $line_num $line_num $i
# 		fi
# 	done
# done


# Fifth round:
log_dir=../logs/rasqual/qsub5/
size=1
echo 5th round
for i in `seq 22`; do
	n_genes=`wc -l ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | cut -d" " -f1`
	for line_num in `seq 1 $size $n_genes`; do
		param=($(cat ../processed_data/rasqual/input/rasqual.input.chr$i.filt.txt | sed "${line_num}q;d"))
		gene_id=${param[0]}
		gene_name=${param[1]}
		out_dir=../processed_data/rasqual/output/galactose/joint/chr$i/
		if [[ -s $out_dir/${gene_id}_${gene_name}.txt ]]; then
			echo "$out_dir/${gene_id}_${gene_name}.txt exist! skipping..."
		else
			echo INFO - chr$i:$line_num 
			qsub -N ras.chr$i.$line_num -l h_rt=99:99:99 -l h_vmem=8G -cwd -o $log_dir/rasqual.chr$i.$line_num.log -e $log_dir/rasqual.chr$i.$line_num.err rasqual/scg4/rasqual.driver.sh $line_num $line_num $i
		fi
	done
done