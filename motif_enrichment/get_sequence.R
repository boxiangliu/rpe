library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)

treeQTL_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
rasqual_dir = c(
	glucose = '../processed_data/rasqual/output/glucose/joint/',
	galactose = '../processed_data/rasqual/output/galactose/joint/'
)
hg19 = '/mnt/lab_data/montgomery/shared/genomes/hg19/hg19.fa'

out_dir = '../processed_data/motif_enrichment/get_sequence/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


read_eGene = function(treeQTL_fn, condition){

	treeQTL = fread(treeQTL_fn)
	setnames(treeQTL, condition, 'target')
	gene_id = treeQTL[target == 1, gene]
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
		lead_SNP_index = which(eQTL$chisq == max_chisq)
		lead_SNP = eQTL[lead_SNP_index,]
		
		return(lead_SNP)

	}

	return(eSNP)

}

get_flanking_sequence = function(chr, pos, reference, window){

	start = pos - window 
	end = pos + window
	command = sprintf('samtools faidx %s %s:%s-%s',reference, chr, start, end)
	x = system(command, intern = TRUE)
	result = toupper(x[2])
	names(result) = x[1]
	return(result)

}

switch_allele = function(sequence, pos, old, new){

	old_len = length(old)
	stopifnot(substr(sequence, pos, pos + old_len - 1) == toupper(old))

	prefix = substr(sequence, 1, pos - 1)
	suffix = substr(sequence, pos + old_len, nchar(sequence))

	new_seq = paste0(prefix,new,suffix)
	names(new_seq) = names(sequence)

	return(new_seq)

}

append_to_fasta = function(sequence, output_fn){

	for (i in seq_along(sequence)){
	
		header = names(sequence[i])
		write(header, output_fn, append=TRUE)
		s = sequence[i]
		write(s,output_fn, append = TRUE)
	
	}

}


for (condition in c('glucose','galactose')){
	# Read eGenes
	eGene = read_eGene(treeQTL_fn, condition)

	# Read lead eQTLs
	dir = rasqual_dir[condition]
	eSNP = read_lead_eQTLs(eGene, dir)

	window = 7
	pos_in_seq = 8
	for (i in 1:nrow(eSNP)){
		print(sprintf('sequence %s',i))
		snp = eSNP[i,]

		chr = snp[,chr]
		pos = snp[,pos]
		
		# get 15-bp around the lead eQTL variant
		ase = snp[,ase]
		ref = snp[,ref]
		alt = snp[,alt]
		if (ase > 0.5){
			background = get_flanking_sequence(chr, pos, reference = hg19, window = window)
			target = switch_allele(background, pos_in_seq, ref, alt)
		} else if (ase < 0.5){
			target = get_flanking_sequence(chr, pos, reference = hg19, window = window)
			background = switch_allele(target, pos_in_seq, ref, alt)
		} else {
			next
		}

		# save the sequences to fasta file
		target_fn = sprintf('%s/%s.target.fa', out_dir, condition)
		background_fn = sprintf('%s/%s.background.fa', out_dir, condition)
		append_to_fasta(target,target_fn)
		append_to_fasta(background,background_fn)
		
	}
}

