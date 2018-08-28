# format_genotype selects the lead eQTL and grabs their 
# genotype from a VCF file. 
library(data.table)
library(doMC)
library(foreach)
registerDoMC(10)
library(stringr)

treeQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
rasqual_dir = c('../processed_data/rasqual/output/glucose/joint/')
vcf_dir = '/srv/persistent/bliu2/rpe/data/genotype/filt/'
out_dir = '../processed_data/phenotype_association/format_genotype/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_eGene = function(treeQTL_fn){

	treeQTL = fread(treeQTL_fn)
	gene_id = treeQTL[, gene]
	return(gene_id)

}

read_lead_eQTLs = function(eGene, dir){

	eSNP = foreach(gene = eGene, .combine = rbind)%dopar%{

		fn = list.files(dir, pattern = gene, full.names = TRUE, recursive = TRUE)
		print(fn)
		stopifnot(length(fn) == 1)

		select = c(1:6,11:12)
		col.names = c('gene','snp','chr','pos','ref','alt','chisq','ase')
		eQTL = fread(fn, select = select, col.names = col.names)
		
		max_chisq = max(eQTL$chisq)
		lead_SNP_index = which(eQTL$chisq == max_chisq)[1] # breaking ties by selecting the first one
		lead_SNP = eQTL[lead_SNP_index,]
		
		return(lead_SNP)

	}

	return(eSNP)

}

read_VCF = function(vcf_dir, chr, pos){

	vcf_fn = list.files(vcf_dir, pattern = paste0(chr,'.all_filters.vcf.gz$'), full.names=TRUE)
	
	chr = str_replace(chr,'chr','')
	command = sprintf("bcftools view %s %s:%s-%s | bcftools query -H -f '%%CHROM\t%%POS\t%%ID[\t%%DS]\n'",vcf_fn, chr, pos, pos)

	result = system(command, intern = TRUE)
	result = str_split(result,pattern='\t')


	header = str_extract(result[[1]][-c(1,2,3)],'(?<=])(.+?)(?=:)')

	print(result[[2]][2])
	print(pos)

	remove_row = c()
	if (length(result) > 2){
		for (i in 2:length(result)){
			if (as.integer(result[[i]][2]) != pos){
				remove_row = c(remove_row, i)
			}
		}
		result = result[-remove_row]
	}
	
	stopifnot(as.integer(result[[2]][2]) == pos)
	
	rsid = result[[2]][3]
	dosage = as.numeric(result[[2]][-c(1,2,3)])

	data = data.table(chr = chr, pos = pos, rsid = rsid, sample = header, dosage = dosage)

	return(data)

}

dna2rna = fread('/srv/persistent/bliu2/rpe/data/meta/dna2rna.txt',colClasses = 'character')
old = unlist(dna2rna[,DNA])
new = unlist(dna2rna[,RNA])

convert_sample_names = function(sample, old, new){

	map = new
	names(map) = old
	return(unname(map[sample]))

}

eGene = read_eGene(treeQTL_fn)
lead_eQTLs = read_lead_eQTLs(eGene, rasqual_dir)

dosage = foreach (i = 1:nrow(lead_eQTLs),.combine = 'rbind')%do%{

	snp = lead_eQTLs[i,]
	chr = snp[,chr]
	pos = snp[,pos]

	read_VCF(vcf_dir, chr, pos)

}

dosage$sample = convert_sample_names(dosage$sample,old,new)

out_fn = sprintf('%s/dosage.txt',out_dir)
fwrite(dosage, out_fn, sep = '\t')