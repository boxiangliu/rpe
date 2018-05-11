ls /srv/persistent/bliu2/rpe/data/rnaseq/fastq/all_bams/ | \
sed 's/_R1.fastq.gz//;s/_R2.fastq.gz//' | \
uniq > rnaseq/samples.txt

while read sample; do
echo $sample
mkdir ../data/rnaseq/bam_gtex_pipe/${sample}/

/srv/persistent/bliu2/tools/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--runThreadN 15 \
--genomeDir /srv/persistent/bliu2/shared/aligner_index/STAR/hg19_gencode19_overhang75 \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--readFilesIn ../data/rnaseq/fastq/all_bams/${sample}_R1.fastq.gz ../data/rnaseq/fastq/all_bams/${sample}_R2.fastq.gz \
--readFilesCommand zcat \
--outFileNamePrefix ../data/rnaseq/bam_gtex_pipe/${sample}/${sample}. \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs None \
--alignSoftClipAtReferenceEnds Yes \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--genomeLoad NoSharedMemory \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimOutType WithinBAM SoftClip \
--outSAMattributes NH HI AS nM NM \
--outSAMattrRGline ID:${sample} SM:${sample}

samtools index ../data/rnaseq/bam_gtex_pipe/${sample}/${sample}.Aligned.sortedByCoord.out.bam

java -Xmx2g -jar /software/picard-tools/1.92/MarkDuplicates.jar \
I=../data/rnaseq/bam_gtex_pipe/${sample}/${sample}.Aligned.sortedByCoord.out.bam \
O=../data/rnaseq/bam_gtex_pipe/${sample}/${sample}.Aligned.sortedByCoord.out.dup.bam \
M=../data/rnaseq/bam_gtex_pipe/${sample}/marked_dup_metrics.txt

samtools index ../data/rnaseq/bam_gtex_pipe/${sample}/${sample}.Aligned.sortedByCoord.out.dup.bam

done < rnaseq/samples.txt
