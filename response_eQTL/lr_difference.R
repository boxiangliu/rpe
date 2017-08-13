library(data.table)
library(stringr)

in_dir='../processed_data/response_eQTL/residual/'

glu=fread(sprintf('%s/glucose_residual.tsv',in_dir))
gal=fread(sprintf('%s/galactose_residual.tsv',in_dir))
stopifnot(dim(glu)==dim(gal))
stopifnot(all(glu$gene_id==gal$gene_id))

diff=glu[,2:ncol(glu)]-gal[,2:ncol(gal)]
setnames(diff,colnames(diff),str_replace(colnames(diff),'\\.glucose',''))
diff$gene_id=glu$gene_id
setcolorder(diff,c(ncol(diff),1:(ncol(diff)-1)))