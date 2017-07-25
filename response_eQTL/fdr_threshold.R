library(data.table)
library(cowplot)

# Variables:
out_dir='../processed_data/response_eQTL/fdr_threshold/'
fig_dir='../figures/response_eQTL/fdr_threshold/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


# Read treeQTL results:
# Glucose:
glu_fn='../processed_data/rasqual/output/glucose/treeQTL/fdr0.5/eAssocs.txt'
glu=fread(glu_fn)
glu[,id:=paste(gene,SNP,sep='_')]

# Galactose:
gal_fn='../processed_data/rasqual/output/galactose/treeQTL/fdr0.5/eAssocs.txt'
gal=fread(gal_fn)
gal[,id:=paste(gene,SNP,sep='_')]

# Find condition-specific eQTL:
# Glucose:
glu_spec=glu[!(id%in%gal$id)&BBFDR<0.001,]
glu_spec[,minFDR:=min(BBFDR),by='gene']
glu_out=glu_spec[BBFDR==minFDR,list(gene,SNP,p.value,BBFDR)]
fwrite(glu_out,sprintf('%s/glucose_specific.txt',out_dir),sep='\t')

# Galactose:
gal_spec=gal[!(id%in%glu$id)&BBFDR<0.001,]
gal_spec[,minFDR:=min(BBFDR),by='gene']
gal_out=gal_spec[BBFDR==minFDR,list(gene,SNP,p.value,BBFDR)]
fwrite(gal_out,sprintf('%s/galactose_specific.txt',out_dir),sep='\t')


# Read RASQUAL result for condition-specific eQTLs:
pattern=paste(unique(c(glu_out$gene,gal_out$gene)),collapse='|')
glu_fn_list=list.files('../processed_data/rasqual/output/glucose/joint/',pattern=pattern,recursive=TRUE,full.names=TRUE)
container=list()
for (f in glu_fn_list){
	temp=fread(f,select=c(1,2,3,4,5,6,11,12,25),col.names=c('gene','SNP','chr','pos','ref','alt','chisq','pi','r2_rsnp'))
	container[[f]]=temp[SNP%in%c(glu_out$SNP,gal_out$SNP)]
}
eqtl_glu_spec=Reduce(rbind,container)


gal_fn_list=list.files('../processed_data/rasqual/output/galactose/joint/',pattern=pattern,recursive=TRUE,full.names=TRUE)
container=list()
for (f in gal_fn_list){
	temp=fread(f,select=c(1,2,3,4,5,6,11,12,25),col.names=c('gene','SNP','chr','pos','ref','alt','chisq','pi','r2_rsnp'))
	container[[f]]=temp[SNP%in%c(glu_out$SNP,gal_out$SNP)]
}
eqtl_gal_spec=Reduce(rbind,container)


stopifnot(all(dim(eqtl_glu_spec)==dim(eqtl_gal_spec)))

eqtl_glu_spec$treatment='glucose'
eqtl_gal_spec$treatment='galactose'

eqtl=rbind(eqtl_glu_spec,eqtl_gal_spec)
eqtl[,pval:=pchisq(chisq,df=1,lower.tail=FALSE)]
eqtl[,logp:=-log10(pval)]
eqtl_cast=dcast(eqtl,gene+SNP~treatment,value.var=c('logp','pi','r2_rsnp'))


# Filter out rSNP with low r2: 
eqtl_cast_filt=eqtl_cast[r2_rsnp_galactose>=0.95&r2_rsnp_glucose>=0.95]
eqtl_cast_filt[,c('logp_diff','pi_diff'):=list(logp_glucose-logp_galactose,pi_glucose-pi_galactose)]


# Plot difference in P-value vs. difference in effect size (pi):
p=ggplot(eqtl_cast_filt,aes(logp_diff,pi_diff,size=abs(logp_diff*pi_diff)))+geom_point(alpha=0.2)+scale_size_continuous(guide='none')+geom_hline(yintercept=c(-0.1,0.1),linetype='dashed',color='red')+geom_vline(xintercept=c(2,-2),linetype='dashed',color='red')+xlab('Difference in -log10(P-value)')+ylab('Difference in effect size')
save_plot(sprintf('%s/diff_p_vs_diff_pi.pdf',fig_dir),p)


# Output result: 
eqtl_cast_filt[,temp:=abs(logp_diff*pi_diff)]
setorder(eqtl_cast_filt,-temp)
eqtl_cast_filt[,c('r2_rsnp_galactose','r2_rsnp_glucose','temp'):=NULL]
fwrite(eqtl_cast_filt,sprintf('%s/condition_specific_eqtl.tsv',out_dir),sep='\t')
