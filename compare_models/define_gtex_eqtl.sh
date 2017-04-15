in_dir=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_allpairs_FOR_QC_ONLY
out_dir=../processed_data/compare_models/gtex_eqtl/
if [[ ! -d $out_dir ]]; then mkdir -p $out_dir; fi

# Subset to chr22: 
for f in $(ls $in_dir/*allpairs.txt.gz); do
	f=$(basename $f)
	tissue=${f/_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz/}
	echo INFO - $tissue: $f
	zgrep -P '\t22_' $in_dir/$f | awk -v tissue=$tissue 'BEGIN{OFS="\t"}{print $0,tissue}' > $out_dir/${f/allpairs.txt.gz/allpairs.chr22.txt}
done


# Concatenate: 
echo INFO - concatenating...
cat $out_dir/*.allpairs.chr22.txt > $out_dir/all_tissues.chr22.txt


# Sort by fid and sid:
echo INFO - sorting...
sort -k1,1 -k2,2 $out_dir/all_tissues.chr22.txt > $out_dir/all_tissues.chr22.sorted.txt


# Select eQTL with lowest p-value:
echo INFO - getting lowest p-value...
cat $out_dir/all_tissues.chr22.sorted.txt | python compare_models/get_lowest_pval.py > $out_dir/all_tissues.chr22.lowest_pval.txt


# Multiple hypothesis correction:
echo INFO - running multiple hypothesis correction with TreeQTL...
Rscript compare_models/treeQTL.R


# Remove intermediate files:
rm $out_dir/*allpairs.chr22.txt
rm $out_dir/all_tissues.chr22.txt
gzip $out_dir/all_tissues.chr22.sorted.txt

