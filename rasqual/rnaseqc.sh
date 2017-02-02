while read line; do
	line=($line)
	in_bam=${line[1]}
	sample=${line[0]}; echo $sample
	java -Xmx2g -jar /software/picard-tools/1.92/AddOrReplaceReadGroups.jar I=$in_bam O=${in_bam/bam/rg.bam} RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> ../logs/AddOrReplaceReadGroups.$sample.log
	samtools index ${in_bam/bam/rg.bam}
done < rasqual/rnaseqc.glucose.txt

while read line; do
	line=($line)
	in_bam=${line[1]}
	sample=${line[0]}; echo $sample
	java -Xmx2g -jar /software/picard-tools/1.92/AddOrReplaceReadGroups.jar I=$in_bam O=${in_bam/bam/rg.bam} RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> ../logs/AddOrReplaceReadGroups.$sample.log
	samtools index ${in_bam/bam/rg.bam}
done < rasqual/rnaseqc.galactose.txt

../../tools/jre1.6.0_45/bin/java -Xmx6g -jar ../../tools/RNA-SeQC_v1.1.8.jar \
-n 1000 -s rasqual/rnaseqc.glucose.txt \
-t ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf \
-r ../../shared/genomes/hg19/hg19.fa \
-o ../data/rnaseq/rnaseqc -noDoC -strictMode > ../logs/rnaseqc.glucose.log


java -Xmx6g -jar ../../tools/RNA-SeQC_v1.1.8.jar \
-n 1000 -s rasqual/rnaseqc.galactose.txt \
-t ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf \
-r ../../shared/genomes/hg19/hg19.fa \
-o ../data/rnaseq/rnaseqc -noDoC -strictMode > ../logs/rnaseqc.galactose.log
