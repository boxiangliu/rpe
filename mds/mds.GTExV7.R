#!/usr/bin/env Rscript
# boxiang liu
# durga
# perform multidimentional scaling 

# library:
library(dplyr)
library(data.table)
library('MASS')
library(cowplot)
library(stringr)
library(ggrepel)

# command line arguments: 
rpkm_file='../processed_data/mds/preprocess.GTExV7/combined.rpkm'
coldata_file='../processed_data/mds/preprocess.GTExV7/combined.col'
color_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gtex_tissue_colors.txt'
figure_dir='../figures/mds/mds.GTExV7/'
out_dir='../processed_data/mds/mds.GTExV7/'
if (!dir.exists(figure_dir)){dir.create(figure_dir)}
if (!dir.exists(out_dir)){dir.create(out_dir)}


# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)
outliers = c('GTEX-Y5LM-0126-SM-4VBRL') # GTEX-Y5LM-0126-SM-4VBRL is an outlier
rpkm[,`GTEX-Y5LM-0126-SM-4VBRL`:=NULL]
col_data = col_data[!(sample %in% outliers)]
stopifnot(all(colnames(rpkm)[2:ncol(rpkm)] == col_data$sample))

# create rpkm matrix:
row_data=rpkm$Name
rpkm_mat=as.matrix(subset(rpkm,select=-Name))
rownames(rpkm_mat)=row_data



# calculate distance based on pearson correlation:
message('calculating distance...')
pearson=cor(rpkm_mat)
dist=as.dist(1-pearson)

# perform multidimensional scaling: 
message('multidimensional scaling...')
mds_res=isoMDS(dist,k=2,maxit=500)
saveRDS(list(pearson,dist,mds_res,rpkm_mat,col_data,row_data),paste0(out_dir,'/mds.rds'))




# find centroid of each tissue: 
mds=data.table(mds_res$points,keep.rownames=T)
colnames(mds)=c('sample','x','y')
mds[x>1,]
mds = mds[x<1,]

mds$tissue=col_data$tissue
# mds[,x_centroid:=median(x),by='tissue']
# mds[,y_centroid:=median(y),by='tissue']
mds[,x_centroid:=mean(x),by='tissue']
mds[,y_centroid:=mean(y),by='tissue']


# prepare color palette: 
color=fread(color_file)[,list(tissue=tissue_site_detail,tissue_color_hex)]
color=rbind(color,
	data.frame(tissue=c('RPE (glu)','RPE (gal)'),
		tissue_color_hex=c('808080','D3D3D3')))
color[,tissue_color_hex:=paste0('#',tissue_color_hex)]


mds=merge(mds,color[,list(color=tissue_color_hex,tissue)],by='tissue')
mds$tissue=as.character(mds$tissue)
mds[,label:=ifelse(tissue%in%c('Liver','Muscle - Skeletal','Brain - Frontal Cortex (BA9)','Whole Blood','Cells - EBV-transformed lymphocytes','Pancreas','Testis','RPE (glu)','Kidney - Cortex'),tissue,'')]
mds[,label:=str_replace(label,' - Frontal Cortex \\(BA9\\)','')]
mds[,label:=str_replace(label,'Cells - EBV-transformed lymphocytes','LCL')]
mds[,label:=str_replace(label,'RPE (glu)','RPE')]
mds[,label:=str_replace(label,'Kidney - Cortex','Kidney')]


# make scatter plot:
message('plotting...')
mds_centroid=unique(mds[,list(x_centroid,y_centroid,tissue,label)])
color_map=color$tissue_color_hex
names(color_map)=color$tissue

p1=ggplot(mds,aes(x=x,y=y,color=tissue))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text(data=mds_centroid,aes(x=x_centroid,
		y=y_centroid,label=label),color='black')+
	theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15))


mds[,label:=ifelse(tissue%in%'RPE (glu)',sample,'')]

p2=ggplot(mds,aes(x=x,y=y,color=tissue))+geom_point()+
	scale_color_manual(values=color_map)+
	geom_text_repel(aes(x=x,y=y,label=label),color='black')+
	theme_bw()+xlab('MDS Coordinate 1')+
	ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15))

pdf(paste0(figure_dir,'/mds.pdf'),height=6,width=14.5)
p1;p2
dev.off()

# Small MDS test:
idx = seq(1,ncol(pearson),by=10)
small_pearson = pearson[idx,idx]
small_dist=as.dist(1-small_pearson)
small_col_data = col_data[match(colnames(small_pearson),col_data$sample),]
stopifnot(all(small_col_data$sample == colnames(small_pearson)))

# perform multidimensional scaling: 
message('multidimensional scaling...')
small_mds_res=isoMDS(small_dist,k=2,maxit=500)
saveRDS(list(small_pearson,small_dist,small_mds_res,small_col_data),paste0(out_dir,'/small_mds.rds'))



# find centroid of each tissue: 
small_mds=data.table(small_mds_res$points,keep.rownames=T)
colnames(small_mds)=c('sample','x','y')
small_mds$tissue=small_col_data$tissue
# mds[,x_centroid:=median(x),by='tissue']
# mds[,y_centroid:=median(y),by='tissue']
small_mds[,x_centroid:=mean(x),by='tissue']
small_mds[,y_centroid:=mean(y),by='tissue']


# prepare color palette: 
color=fread(color_file)[,list(tissue=tissue_site_detail,tissue_color_hex)]
color=rbind(color,
	data.frame(tissue=c('RPE (glu)','RPE (gal)'),
		tissue_color_hex=c('808080','D3D3D3')))
color[,tissue_color_hex:=paste0('#',tissue_color_hex)]


small_mds=merge(small_mds,color[,list(color=tissue_color_hex,tissue)],by='tissue')
small_mds$tissue=as.character(small_mds$tissue)

small_mds[,label:=ifelse(tissue%in%c('Liver','Muscle - Skeletal','Brain - Frontal Cortex (BA9)','Whole Blood','Cells - EBV-transformed lymphocytes','Pancreas','Testis','RPE (gal)','Kidney - Cortex','Heart - Atrial Appendage','Heart - Left Ventricle','Cells - Transformed fibroblasts','Lung','Pituitary','Skin - Not Sun Exposed (Suprapubic)','Esophagus - Mucosa'),tissue,'')]
small_mds[,label:=str_replace(label,' - Frontal Cortex \\(BA9\\)','')]
small_mds[,label:=str_replace(label,'Cells - EBV-transformed lymphocytes','LCL')]
small_mds[,label:=str_replace(label,'RPE (gal)','RPE')]
small_mds[,label:=str_replace(label,'Kidney - Cortex','Kidney')]
small_mds[,label:=str_replace(label,'Heart - Atrial Appendage','Heart')]
small_mds[,label:=str_replace(label,'Heart - Left Ventricle','Heart')]
small_mds[,label:=str_replace(label,'Cells - Transformed fibroblasts','Fibroblasts')]
small_mds[,label:=str_replace(label,'Skin - Not Sun Exposed \\(Suprapubic\\)','Skin')]
small_mds[,label:=str_replace(label,'Esophagus - Mucosa','Esophagus')]


# make scatter plot:
message('plotting...')
small_mds_centroid=unique(small_mds[,list(x_centroid,y_centroid,tissue,label)])
color_map=color$tissue_color_hex
names(color_map)=color$tissue

p1=ggplot(small_mds,aes(x=x,y=y,color=tissue))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text(data=small_mds_centroid,aes(x=x_centroid,
		y=y_centroid,label=label),color='black')+
	theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15))


p2=ggplot(small_mds,aes(x=x,y=y,color=tissue))+geom_point()+
	scale_color_manual(values=color_map)+
	theme_bw()+xlab('MDS Coordinate 1')+
	ylab('MDS Coordinate 2')+
	theme(axis.title=element_text(size=15))

pdf(paste0(figure_dir,'/small_mds.pdf'),height=6,width=14.5)
p1;p2
dev.off()

# Calculate PCs:
pc_res=prcomp(t(rpkm_mat))
saveRDS(list(pc_res,rpkm_mat,col_data,row_data),paste0(out_dir,'/pc.rds'))
pc=data.table(pc_res$x[,1:10],keep.rownames=T)
colnames(pc)=c('sample',paste0('PC',seq(ncol(pc)-1)))
pc$tissue=col_data$tissue


for (i in seq(ncol(pc)-2)){
	col1=paste0('PC',i)
	col2=paste0('PC',i,'_centroid')
	pc[,c(col2):=mean(get(col1)),by='tissue']
}


# prepare color palette: 
color=fread(color_file)[,list(tissue=tissue_site_detail,tissue_color_hex)]
color=rbind(color,
	data.frame(tissue=c('RPE (glu)','RPE (gal)'),
		tissue_color_hex=c('808080','D3D3D3')))
color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
pc=merge(pc,color[,list(tissue,color=tissue_color_hex)],by='tissue')
pc$tissue=as.character(pc$tissue)

# make scatter plot:
message('plotting...')
color_map=color$tissue_color_hex
names(color_map)=color$tissue

pc[,label:=ifelse(tissue%in%'RPE (glu)',sample,'')]

p4=ggplot(pc,aes(x=PC1,y=PC2,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p5=ggplot(pc,aes(x=PC1,y=PC3,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p6=ggplot(pc,aes(x=PC2,y=PC3,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p7=ggplot(pc,aes(x=PC3,y=PC4,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))


p8=ggplot(pc,aes(x=PC4,y=PC5,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p9=ggplot(pc,aes(x=PC5,y=PC6,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p10=ggplot(pc,aes(x=PC6,y=PC7,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p11=ggplot(pc,aes(x=PC7,y=PC8,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))


p12=ggplot(pc,aes(x=PC8,y=PC9,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))

p13=ggplot(pc,aes(x=PC9,y=PC10,color=tissue,label=label))+
	geom_point()+scale_color_manual(values=color_map)+
	geom_text_repel(color='black')+
	theme_bw()+
	theme(axis.title=element_text(size=15))


pdf(paste0(figure_dir,'/pc.pdf'),width=14.5,height=6)
p4;p5;p6;p7;p8;p9;p10;p11;p12;p13
dev.off()
