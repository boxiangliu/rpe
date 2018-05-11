rnaseq/collapse_annotation.py \
--transcript_blacklist rnaseq/gencode19_unannotated_readthrough_blacklist.txt \
/mnt/lab_data/montgomery/shared/annotations/gencode.v19.annotation.gtf \
../data/reference/gencode.v19.annotation.collapsed_annotation.gtf

run_rnaseqc(){

sample=$1
echo $sample
mkdir -p ../data/rnaseq/rnaseqc/${sample}/

/srv/persistent/bliu2/tools/jre1.7.0_80/bin/java -Xmx6g \
-jar /srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar \
-n 1000 \
-s "${sample}|../data/rnaseq/bam_gtex_pipe/${sample}/${sample}.Aligned.sortedByCoord.out.dup.bam|${sample}" \
-t ../data/reference/gencode.v19.annotation.collapsed_annotation.gtf \
-r /srv/persistent/bliu2/shared/genomes/hg19/hg19.fa \
-o ../data/rnaseq/rnaseqc/${sample}/ \
-noDoC -strictMode 

}

export -f run_rnaseqc
parallel -j15 run_rnaseqc {} :::: rnaseq/samples.txt
