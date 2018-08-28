################
# Get sequence #
################
This analysis will test if the lead variants in eQTLs from each condition are enriched in any motifs, and to see if the enriched motifs are bound by TFs differentially expressed across two conditions.

########
# Main #
########

for condition in {glucose,galactose}:
	
	Read eGenes
	
	Read lead eQTLs

	for SNP in lead eQTLs:

		sequence = get 15-bp around the lead eQTL variant;

		if direction is up-regulation:

			target sequence = change the lead eQTL variant to alt allele;

			background sequence = sequence;

		else:

			target sequence = sequence;

			background sequence = change the lead eQTL variant to alt allele;

		end if

	end for 

	save the target sequences to fasta file 

	save the background sequences to fasta file 
	
end for

###########
# Methods #
###########
declare treeQTL_fn
declare rasqual_dir

method: read eGenes(treeQTL_fn, condition):
	
	read treeQTL_fn

	subset to condition

	return eGenes

end method


method: read lead eQTL(eGenes, rasqual_dir):

	for gene in eGenes:

		get lead eSNP for eGenes 

	end for 

end method


method: get lead eSNP for gene(rasqual_dir, gene):
	
	read the eQTL file for the gene from rasqual_dir

	get the lead eSNP 

	return {gene, SNP, chrom, position, ref, alt, direction}

end method

declare reference = the reference genome

method: get flanking sequence(SNP, reference, window):

	get flanking sequence of window size 

end method


method: switch allele(sequence, position, old, new):

	change the letter in sequence:position from old to new

end method


declare output_fn 

method: save_fasta(list of sequences, output_fn):
	
	for sequence in the list of sequences: 

		format string 

		write string to file 

	end for 

end method 


declare target_fn

declare background_fn

declare output_dir

######### 
# HOMER #
#########
method: perform HOMER analysis(target_fn, background_fn, output_dir):
	
	HOMER target_fn output_dir background_fn

end method  
