library(data.table)
library(foreach)
library(doMC)
registerDoMC(15)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(openxlsx)


eGenes_fn = '../processed_data/response_eQTL/treeQTL_MT/Bliu_MTtreeQTL/eGenesMT.txt'
glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint'
genotype_pc_fn = '../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv'
eyegex_fn = '../data/eyegex/eyegex-supp3.xlsx'

out_dir = '../processed_data/genotype_pc_eQTL_correlation/correlation/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive = TRUE)}

fig_dir = '../figures/genotype_pc_eQTL_correlation/correlation/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive = TRUE)}

# Functions
read_lead_eQTLs = function(eGene, dir){

	eSNP = foreach(gene = eGene, .combine = rbind)%dopar%{

		fn = list.files(dir, pattern = gene, full.names = TRUE, recursive = TRUE)
		print(fn)
		stopifnot(length(fn) == 1)

		select = c(1:6,11:12)
		col.names = c('gene','snp','chr','pos','ref','alt','chisq','ase')
		eQTL = fread(fn, select = select, col.names = col.names)
		
		max_chisq = max(eQTL$chisq)
		lead_SNP_index = sample(which(eQTL$chisq == max_chisq),1) # Break ties at random
		lead_SNP = eQTL[lead_SNP_index,]
		return(lead_SNP)

	}
	eSNP$chr = str_replace(eSNP$chr, 'chr', '')
	res = str_split_fixed(eSNP$gene,'_',2)
	eSNP$gene = res[,1]
	eSNP$gene_name = res[,2]
	eSNP$gene = str_split_fixed(eSNP$gene,'\\.',2)[,1]

	return(eSNP)
}


extract_lead_eSNP_genotype = function(eSNP){

	region_fn = paste0(out_dir,'lead_SNPs.txt')
	fwrite(eSNP[,list(chr,pos,pos)],region_fn,sep='\t',col.names = FALSE)
	out_fn = paste0(out_dir,'genotypes.txt')

	for (i in seq(22)){
		print(i)
		vcf_fn = sprintf('../data/genotype/filt/rpe.imputed.chr%s.all_filters.vcf.gz',i)

		if (i == 1){
			cmd = sprintf('bcftools view -R %s %s | bcftools query -H -f "%%CHROM\\t%%POS[\\t%%DS]\\n" > %s', region_fn, vcf_fn, out_fn)
		} else {
			cmd = sprintf('bcftools view -R %s %s | bcftools query -f "%%CHROM\\t%%POS[\\t%%DS]\\n" >> %s', region_fn, vcf_fn, out_fn)
		}
		system(cmd)
	}

	genotypes = fread(out_fn)
	colnames(genotypes) = paste0(colnames(genotypes),':')
	colnames(genotypes) = str_extract(colnames(genotypes),'(?<=])(.+?)(?=:)')
	genotypes$CHROM = as.character(genotypes$CHROM)

	unlink(region_fn)
	unlink(out_fn)

	return(genotypes)
}

# Get RPE eGenes 
eGenes = fread(eGenes_fn)
eGenes$gene = str_split_fixed(eGenes$gene,'\\.',2)[,1]

# Get RPE lead eSNP
glucose_eSNP = read_lead_eQTLs(eGenes[glucose==1,gene],glucose_dir)
galactose_eSNP = read_lead_eQTLs(eGenes[galactose==1,gene],galactose_dir)

# Get lead eSNP genotype 
glucose_genotype = extract_lead_eSNP_genotype(glucose_eSNP)
galactose_genotype = extract_lead_eSNP_genotype(galactose_eSNP)


# Get genotype PCs 
genotype_pc = fread(genotype_pc_fn)
dna2rna = fread('../data/meta/dna2rna.txt', colClasses = 'character')
genotype_pc$sample = dna2rna[match(genotype_pc$sample, dna2rna$RNA),DNA]


# Read eyegex dataset: 
eyegex = as.data.table(read.xlsx(eyegex_fn,sheet = 1, rows = 6:14865))



calc_max_r2 = function(genotype,genotype_pc,eSNP,eyegex){
	# Calculate correlation between genotype PCs and eSNP genotype
	mat1 = genotype[,match(genotype_pc$sample,colnames(genotype)),with=FALSE]
	mat2 = genotype_pc[,2:11]


	# Calculate the mean of each vector
	mat1_means = rowMeans(mat1)
	mat2_means = colMeans(mat2)

	# Subtract the mean from each vector 
	mat1 = sweep(mat1,1,mat1_means)
	mat2 = sweep(mat2,2,mat2_means)

	# Calculate the standard deviation of each vector 
	mat1_sd = apply(mat1,1,sd)
	mat2_sd = apply(mat2,2,sd)

	# Divide each vector by its S.D. 
	mat1 = sweep(mat1,1,mat1_sd,FUN = '/')
	mat2 = sweep(mat2,2,mat2_sd,FUN = '/')

	# Multiply standardized genotype and PC matrices and take squares
	r2 = (as.matrix(mat1) %*% as.matrix(mat2) / 23)^2


	# Take the maximum of squared matrix correlation coefficients
	r2_max = apply(r2,1,max)
	r2_max = data.table(genotype[,list(chr = CHROM, pos = POS)], r2_max = r2_max)
	r2_max = merge(r2_max,eSNP[,list(gene,gene_name,chr,pos)],by=c('chr','pos'))
	r2_max = merge(r2_max,eGenes,by=c('gene','gene_name'))
	r2_max$eyegex = r2_max$gene %in% eyegex$gene

	return(r2_max)
}


glucose_r2_max = calc_max_r2(glucose_genotype,genotype_pc,glucose_eSNP,eyegex)
# x1 = glucose_r2_max[,list(r2_max, condition = 'glucose', gene, specific = ifelse(galactose == 1, 'Shared','Specific'))] # condition-specific eQTL
x1 = glucose_r2_max[,list(r2_max, condition = 'glucose', gene, specific = ifelse(eyegex == TRUE, 'Shared with EyeGEx','RPE-specific'))] 

galactose_r2_max = calc_max_r2(galactose_genotype,genotype_pc,galactose_eSNP,eyegex)
# x2 = galactose_r2_max[,list(r2_max, condition = 'galactose', gene, specific = ifelse(glucose == 1, 'Shared', 'Specific'))]
x2 = galactose_r2_max[,list(r2_max, condition = 'galactose', gene, specific = ifelse(eyegex == TRUE, 'Shared with EyeGEx','RPE-specific'))]


merged_r2_max = rbind(x1,x2)


# Stratify eGenes into shared and condition-specific
p = ggplot(merged_r2_max,aes(x = as.character(specific), y = r2_max)) + 
	geom_boxplot() + 
	xlab('eQTL classification') + 
	ylab(expression('max('~r^2~')')) + 
	geom_signif(comparisons = list(c('Shared with EyeGEx','RPE-specific')))


t.test(r2_max~specific,merged_r2_max) # p = 0.7334
wilcox.test(r2_max~specific,merged_r2_max) # p = 0.6937

fig_fn = sprintf('%s/max_r2_comparision.pdf', fig_dir)
save_plot(fig_fn, p)
