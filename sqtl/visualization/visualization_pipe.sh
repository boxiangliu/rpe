export leafcutter='/srv/persistent/bliu2/tools/leafcutter/'
gtf2leafcutter_dir='../processed_data/sqtl/visualization/gtf2leafcutter/' 
mkdir -p $gtf2leafcutter_dir
$leafcutter/leafviz/gtf2leafcutter.pl -o $gtf2leafcutter_dir ../data/reference/gencode.v19.annotation.collapsed_annotation.gtf

mkdir -p ../processed_data/sqtl/visualization/prepare_results/

zcat ../data/rnaseq/leafcutter/glucose/cluster/sqtl_perind_numers.counts.gz | \
grep -e "chr12:56115278:56115473:clu_1202" \
-e "chr12:56115278:56117670:clu_1202" \
-e "chr12:56115731:56117670:clu_1202" \
-e "021010" | \
gzip > ../processed_data/sqtl/visualization/prepare_results/rdh5_perind_numers.counts.gz

$leafcutter/leafviz/prepare_results.R \
../processed_data/sqtl/visualization/prepare_results/rdh5_perind_numers.counts.gz \
../processed_data/sqtl/visualization/prepare_input/significance.txt \
../processed_data/sqtl/visualization/prepare_input/effect_size.txt \
$gtf2leafcutter_dir \
--output=../processed_data/sqtl/visualization/prepare_results/sqtl.RData \
--meta_data_file=../processed_data/sqtl/visualization/prepare_input/meta_data.txt \
--FDR=0.05

cd $leafcutter/leafviz/
./run_leafviz.R /srv/persistent/bliu2/rpe/processed_data/sqtl/visualization/prepare_results/sqtl.RData