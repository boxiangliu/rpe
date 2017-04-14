library(data.table)

in_fn='/srv/persistent/bliu2/shared/1000genomes/phase3v5a/integrated_call_samples_v3.20130502.ALL.panel'
out_dir='../processed_data/genotype_pc/1kg_samples/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=T)}

set.seed(42)
samples=panel[,list(sample=sample(sample,4,replace=FALSE)),by='pop']
samples=merge(samples,panel,by=c('sample','pop'),sort=FALSE)

fwrite(samples,sprintf('%s/4_sample_each_pop.panel',out_dir),sep='\t')