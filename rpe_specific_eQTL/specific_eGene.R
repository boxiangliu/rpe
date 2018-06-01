library(data.table)
library(VennDiagram)
library(grid)
library(cowplot)

out_dir='../processed_data/rpe_specific_eQTL/specific_eGene/'
fig_dir='../figures/rpe_specific_eQTL/specific_eGene/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

glu_egene_fn='../processed_data/rasqual/output/glucose/treeQTL/eGenes.txt'
glu_egene=fread(glu_egene_fn)

gal_egene_fn='../processed_data/rasqual/output/galactose/treeQTL/eGenes.txt'
gal_egene=fread(gal_egene_fn)

gtex_egene_fn='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p.egenes.auto.mlinc.bh.txt'
gtex_egene=fread(gtex_egene_fn)

annotation_fn='/srv/persistent/bliu2/shared/annotation/gtex/gencode.v19.genes.v6p.hg19.bed'
annotation=fread(annotation_fn,select=c(5,6,7),col.names=c('gene_id','gene_name','type'))
id2name=annotation[type%in%c('protein_coding','lincRNA'),gene_name]
names(id2name)=annotation[type%in%c('protein_coding','lincRNA'),gene_id]


glu_spec_egene_idx=which(!(glu_egene[,family] %in% unique(gtex_egene[fdr<0.05,gene_id])))
glu_spec_egene=glu_egene[glu_spec_egene_idx,]
glu_spec_egene[,gene_name:=id2name[family]]


gal_spec_egene_idx=which(!(gal_egene[,family] %in% unique(gtex_egene[fdr<0.05,gene_id])))
gal_spec_egene=gal_egene[gal_spec_egene_idx,]
gal_spec_egene[,gene_name:=id2name[family]]


ind_mat=data.table(gene_name=id2name,gene_id=names(id2name))
ind_mat[,c('gtex','glu','gal'):=list(
	gene_id%in%gtex_egene[fdr<0.05,gene_id],
	gene_id%in%glu_egene$family,
	gene_id%in%gal_egene$family)]



pdf(sprintf('%s/specific_eGene.pdf',fig_dir))
draw.triple.venn(area1=sum(ind_mat[,glu]),
	area2=sum(ind_mat[,gal]),
	area3=sum(ind_mat[,gtex]),
	n12=sum(ind_mat[,glu&gal]),
	n23=sum(ind_mat[,gal&gtex]),
	n13=sum(ind_mat[,glu&gtex]),
	n123=sum(ind_mat[,glu&gal&gtex]),
	category = c('Glucose','Galactose','GTEx'),
	lty=rep('blank',3),
	fill=c('#DDE2F4','#D1E0D4','orchid'),
	cex=rep(1.5, 7),
	alpha=rep(0.8, 3),
	fontfamily=rep("sans",7),
	cat.fontfamily = rep("sans", 3),
	cat.cex=rep(1.5, 3))
dev.off()

glu_spec_egene=ind_mat[glu&!gal&!gtex,list(gene_id,gene_name)]
glu_spec_egene=merge(glu_egene,glu_spec_egene,by.x='family',by.y='gene_id')
setorder(glu_spec_egene,fam_p)
fwrite(glu_spec_egene,sprintf('%s/glu_spec_egene.txt',out_dir),sep='\t')

gal_spec_egene=ind_mat[!glu&gal&!gtex,list(gene_id,gene_name)]
gal_spec_egene=merge(gal_egene,gal_spec_egene,by.x='family',by.y='gene_id')
setorder(gal_spec_egene,fam_p)
fwrite(gal_spec_egene,sprintf('%s/gal_spec_egene.txt',out_dir),sep='\t')


rpe_egene=ind_mat[glu&gal&!gtex,list(gene_id,gene_name)]
rpe_egene=merge(glu_egene,rpe_egene,by.x='family',by.y='gene_id')
rpe_egene=merge(gal_egene,rpe_egene,by='family',suffix=c('_glu','_gal'))
setorder(rpe_egene,fam_p_glu,fam_p_gal)
fwrite(rpe_egene,sprintf('%s/rpe_egene.txt',out_dir),sep='\t')


shared_egene=ind_mat[gtex&(glu|gal),list(gene_id,gene_name)]
shared_egene=merge(glu_egene,shared_egene,by.x='family',by.y='gene_id',all.y=TRUE)
shared_egene=merge(gal_egene,shared_egene,by='family',suffix=c('_glu','_gal'),all.y=TRUE)
setorder(shared_egene,fam_p_glu,fam_p_gal)

data=rbind(data.table(pval=c(shared_egene$fam_p_glu,shared_egene$fam_p_gal),type='shared'),
	data.table(pval=glu_spec_egene$fam_p,type='glucose'),
	data.table(pval=gal_spec_egene$fam_p,type='galactose'),
	data.table(pval=c(rpe_egene$fam_p_gal,rpe_egene$fam_p_glu),type='rpe'))
data=data[!is.na(pval),]
data[,logp:=-log10(pval)]
p1=ggplot(data,aes(x=type,y=logp))+geom_boxplot()+scale_y_log10()+xlab('eQTL Category')+ylab('-log10(P-value)')
save_plot(sprintf('%s/pval_distribution.pdf',fig_dir),p1)
