#!/bin/bash

# paths:
wd=$1
# wd=/srv/persistent/bliu2/rpe/data/rnaseq/bam/glucose
cd $wd 

out_dir=$2
# out_dir=/srv/persistent/bliu2/rpe/data/rnaseq/rnaseqc/glucose
mkdir -p $out_dir

# make sample file for RNA-seQC:
ls -d *.bam > sample_file.tmp
cat sample_file.tmp | awk -v wd=$wd 'BEGIN {OFS="\t"} {fn=$1;gsub(/_dedup.rg.bam/,"",$1); print $1,wd"/"fn,"glucose"}' > sample_file.txt
rm sample_file.tmp


# run RNA-seQC: 
parallel -j1 bash /srv/persistent/bliu2/rpe/scripts/rnaseq/rnaseqc.core.sh {} $out_dir :::: sample_file.txt

