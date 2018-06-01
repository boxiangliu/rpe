#!/bin/bash
# Nathan Abell
# Script to execute GSEA2 Pre-Ranked on a RNK-format file
# Usage: bash GSEAPreRanked.sh <RNK_file> <output_directory>


# Variable declarations:
export GSEA="/users/bliu2/tools/gsea/"
export out_dir='../processed_data/rpe_specific_genes/GSEA/'
mkdir -p $out_dir

# Fucntions:
GSEAPreranked(){

RNK=$1
OUTDIR=$2
LABEL=$3
GMX=$4
NPERM=$5
LOG=$OUTDIR/GSEA_run.log

java -Xmx32G -cp $GSEA/gsea-3.0.jar \
xtools.gsea.GseaPreranked \
-rnd_seed 2017 \
-rnk $RNK \
-gmx $GMX \
-out $OUTDIR \
-rpt_label  $LABEL \
-nperm $NPERM \
-plot_top_x 50 \
-collapse false > $LOG 2>&1

}

export -f GSEAPreranked


# Main: 
cat ../processed_data/rpe_specific_genes/rpe_specific_genes.GTExV7/all_genes.RPE.txt | awk 'BEGIN{OFS="\t"}{print $10,$5}' | grep -v zscore | grep -v chr > $out_dir/all_genes.RPE.rnk

GSEAPreranked \
$out_dir/all_genes.RPE.rnk \
$out_dir \
RPE_hallmark_biocarta_kegg_reactome_GO \
$GSEA/MSigDB/hallmark_biocarta_kegg_reactome_GO.v6.1.symbols.gmt \
1000

# Select GO terms with q-value = 0 to rerun:
Rscript rpe_specific_genes/select_GO_terms.R
grep -f ../processed_data/rpe_specific_genes/select_GO_terms/rerun.txt \
$GSEA/MSigDB/hallmark_biocarta_kegg_reactome_GO.v6.1.symbols.gmt > \
../processed_data/rpe_specific_genes/select_GO_terms/rerun.gmt

# Rerun GO terms with q-value = 0:
GSEAPreranked \
$out_dir/all_genes.RPE.rnk \
$out_dir \
RPE_GO_rerun \
../processed_data/rpe_specific_genes/select_GO_terms/rerun.gmt \
100000

# Run GSEA using GO terms with 10,000 permutation:
GSEAPreranked \
$out_dir/all_genes.RPE.rnk \
$out_dir \
RPE_GO \
$GSEA/MSigDB/c5.all.v6.1.symbols.gmt \
10000

# Make table Name (NES,FWER p-value): 
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_fibroblast.GseaPreranked.*/gsea_report_for_na_pos_*.xls
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_artery.GseaPreranked.*/gsea_report_for_na_pos_*.xls
# awk '{FS="\t";OFS=""}{printf "%s (%.02f,%.02f)\n",$1,$6,$8}' ../processed_data/differential_expression/GSEA/meta_heart.GseaPreranked.*/gsea_report_for_na_pos_*.xls
