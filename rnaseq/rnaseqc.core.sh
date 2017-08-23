rnaseqc=/srv/persistent/bliu2/tools/RNA-SeQC_v1.1.8.jar
# gencode19=/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf
gencode19=/srv/persistent/bliu2/shared/shared/annotations/gencode.v19.annotation.gtf 
hg19=/srv/persistent/bliu2/shared/genomes/ucsc_hg19/genome.fa
# rRNA=/srv/persistent/bliu2/shared/genomes/rRNA/human_all_rRNA.fasta

line=($1)
out_dir=$2

sample=${line[0]}
bam=${line[1]}
note=${line[2]}

echo INFO - sample: $sample 
echo INFO - bam: $bam 
echo INFO - note: $note

mkdir -p $out_dir/$sample/

/srv/persistent/bliu2/tools/jre1.6.0_45/bin/java -Xmx6g -jar $rnaseqc \
	-n 1000 -s "$sample|$bam|$note" \
	-t $gencode19 -r $hg19 \
	-o $out_dir/$sample/ \
	-strictMode > /srv/persistent/bliu2/rpe/logs/rnaseqc/rnaseqc.$sample.log