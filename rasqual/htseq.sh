# ../../tools/jre1.6.0_45/bin/java -Xmx6g -jar ../../tools/RNA-SeQC_v1.1.8.jar -n 1000 -s rasqual/rnaseqc.glucose.txt -t ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf -r ../../shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -o ../data/rnaseq/rnaseqc/glucose -noDoC -strictMode > ../logs/rnaseqc.glucose.log
# ../../tools/jre1.6.0_45/bin/java -Xmx6g -jar ../../tools/RNA-SeQC_v1.1.8.jar -n 1000 -s rasqual/rnaseqc.galactose.txt -t ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf -r ../../shared/ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -o ../data/rnaseq/rnaseqc/galactose -noDoC -strictMode > ../logs/rnaseqc.galactose.log

# Glucose:  
cut -f1 rasqual/htseq.glucose.txt > rasqual/htseq.glucose.sample.txt
cut -f2 rasqual/htseq.glucose.txt > rasqual/htseq.glucose.bam.txt
parallel -j10 --xapply htseq-count -f bam -s reverse -r pos -a 10 -t exon -i gene_id {1} ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf '>' ../data/rnaseq/count/glucose/{2}.glucose.gene_count :::: rasqual/htseq.glucose.bam.txt :::: rasqual/htseq.glucose.sample.txt
parallel -j10 --xapply htseq-count -f bam -s reverse -r pos -a 10 -t exon -i exon_id {1} ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf '>' ../data/rnaseq/count/glucose/{2}.glucose.exon_count :::: rasqual/htseq.glucose.bam.txt :::: rasqual/htseq.glucose.sample.txt
rm rasqual/htseq.glucose.sample.txt rasqual/htseq.glucose.bam.txt

# Galactose: 
cut -f1 rasqual/htseq.galactose.txt > rasqual/htseq.galactose.sample.txt
cut -f2 rasqual/htseq.galactose.txt > rasqual/htseq.galactose.bam.txt
parallel -j10 --xapply htseq-count -f bam -s reverse -r pos -a 10 -t exon -i gene_id {1} ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf '>' ../data/rnaseq/count/galactose/{2}.galactose.gene_count :::: rasqual/htseq.galactose.bam.txt :::: rasqual/htseq.galactose.sample.txt
parallel -j10 --xapply htseq-count -f bam -s reverse -r pos -a 10 -t exon -i exon_id {1} ../../shared/annotation/gtex/gencode.v19.genes.v6p.hg19.gtf '>' ../data/rnaseq/count/galactose/{2}.galactose.exon_count :::: rasqual/htseq.galactose.bam.txt :::: rasqual/htseq.galactose.sample.txt
rm rasqual/htseq.galactose.sample.txt rasqual/htseq.galactose.bam.txt


# Merge: 
Rscript rasqual/htseq.merge.R 