library(stringr)

extract_dosage = function(region,vcf){
	command = sprintf('bcftools view -r %s %s | bcftools query -H -f "%%CHROM\t%%POS\t%%REF\t%%ALT[\t%%DS]\n"',region,vcf)
	results = system(command,intern=TRUE)
	for (i in 1:length(results)){
		if (i == 1){
			header = str_split(results[i],'\t')[[1]]
			header = paste0(header,':')
			header = str_extract(header,'(?<=])(.+?)(?=:)')
			dosage = data.frame()
		} else {
			x = str_split(results[i],'\t')[[1]]
			dosage = rbind(dosage,x)
		}
	}
	colnames(dosage)=header
	return(dosage)
}