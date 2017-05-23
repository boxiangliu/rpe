addreadgroup(){
	in_bam=$1
	sample=$2
	echo $sample
	if [[ ! -e ${in_bam/bam/rg.bam} ]]; then 
		java -Xmx2g -jar /software/picard-tools/1.92/AddOrReplaceReadGroups.jar I=$in_bam O=${in_bam/bam/rg.bam} RGID=$sample RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$sample 2> ../logs/AddOrReplaceReadGroups.$sample.log
		samtools index ${in_bam/bam/rg.bam}
	else
		echo ${in_bam/bam/rg.bam} exists!
	fi 
}


while read line; do
	line=($line)
	in_bam=${line[1]}
	sample=${line[0]}
	addreadgroup $in_bam $sample
done < rasqual/addreadgroup.glucose.txt

while read line; do
	line=($line)
	in_bam=${line[1]}
	sample=${line[0]}
	addreadgroup $in_bam $sample
done < rasqual/addreadgroup.galactose.txt
