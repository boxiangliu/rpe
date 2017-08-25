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
rpkm_file='../processed_data/mds/preprocess.all_tissue/combined.rpkm'
coldata_file='../processed_data/mds/preprocess.all_tissue/combined.col'
color_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gtex_tissue_colors.txt'
figure_dir='../figures/mds/mds.all_tissue/'
out_dir='../processed_data/mds/mds.all_tissue/'
if (!dir.exists(figure_dir)){dir.create(figure_dir)}
if (!dir.exists(out_dir)){dir.create(out_dir)}


# read input: 
rpkm=fread(rpkm_file,header=T)
col_data=fread(coldata_file,header=T)


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
mds_res=isoMDS(dist,k=2)
saveRDS(list(pearson,dist,mds_res,rpkm_mat,col_data,row_data),paste0(out_dir,'/mds.rds'))



# find centroid of each tissue: 
mds=data.table(mds_res$points,keep.rownames=T)
colnames(mds)=c('sample','x','y')
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
