library(data.table)
library(stringr)

# Variables
out_dir='../processed_data/select_covariates/merge_covariates'
if (!dir.exists(out_dir)){dir.create(out_dir,recursive=TRUE)}

args=commandArgs(T)
splice_pc_fn=args[1]

sex_fn='../processed_data/sex/sex/gender.tsv'
geno_pc_fn='../processed_data/genotype_pc/genotype_pc/pc_nodup.tsv'
dna2rna_fn='../data/meta/dna2rna.txt'

# glu_sva_fn='../processed_data/hidden_covariates/sva/glucose_LibSizCorr_SVAFact4_MeanGeq10_ZeroLeq20.txt'
# gal_sva_fn='../processed_data/hidden_covariates/sva/galactose_LibSizCorr_SVAFact5_MeanGeq10_ZeroLeq20.txt'


# Read data:
sex=fread(sex_fn,colClasses=c('character','character'))
geno_pc=fread(geno_pc_fn)
glu_sva=fread(glu_sva_fn)
gal_sva=fread(gal_sva_fn)
dna2rna=fread(dna2rna_fn,colClasses=c('character','character'))


# Reformatting:
sex=merge(sex,dna2rna,by.x='sample',by.y='DNA')
sex[,c('sample','RNA'):=list(RNA,NULL)]
sex[,gender:=ifelse(gender=='F',0,1)]

setnames(glu_sva,c('sample',paste0('sva',1:4)))
glu_sva[,sample:=str_replace_all(sample,'X|\\.glucose','')]
glu_sva[,sample:=str_replace_all(sample,'\\.','-')]

setnames(gal_sva,c('sample',paste0('sva',1:5)))
gal_sva[,sample:=str_replace_all(sample,'X|\\.galactose','')]
gal_sva[,sample:=str_replace_all(sample,'\\.','-')]


# Merge:
for (condition in c('glu','gal')){
	for (method in c('sva')){
		x=get(sprintf('%s_%s',condition,method))
		merged=merge(merge(sex,geno_pc[,list(sample,PC1,PC2,PC3)],by='sample'),x,by='sample')
		fwrite(merged,sprintf('%s/%s_sex_geno_%s.tsv',out_dir,condition,method),sep='\t')
	}
}

