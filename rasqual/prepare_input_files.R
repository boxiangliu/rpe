library(data.table)
library(dplyr)
library(dtplyr)
library(rasqualTools)

# Read HTSeq count matrix: 
count=fread('../data/rnaseq/count/rpe_htseqcount_matrix.txt')
count_df=as.data.frame(count[,2:ncol(count)])
rownames(count_df)=unlist(count[,1])


# Split count matrix into glucose and galactose: 
glu=count_df%>%select(contains('lucose'))
gla=count_df%>%select(contains('lactose'))


# Save count matrix in both binary format and plain text format: 
saveRasqualMatrices(list(glucose=glu),"../processed_data/rasqual/expression/", file_suffix = "expression")
saveRasqualMatrices(list(galactose=gla),"../processed_data/rasqual/expression/", file_suffix = "expression")


# Read gene GC content: 
gc=fread('../processed_data/rasqual/gcc.exon.txt',header=F)
setnames(gc,c('gene_id','percentage_gc_content'))


size_factors_glu=rasqualCalculateSampleOffsets(glu,gc,gc_correct = TRUE)
sort(gc[duplicated(gc$gene_id),gene_id])