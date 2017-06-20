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
rpkm_file='../processed_data/mds/preprocess/combined.rpkm'
coldata_file='../processed_data/mds/preprocess/combined.col'
color_file='/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2015-01-12/reference_files/gtex_tissue_colors.txt'
figure_dir='../figures/mds/mds/'
out_dir='../processed_data/mds/mds/'
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
color=fread(color_file)[,list(`Tissues`=tissue_site_detail,tissue_site_detail_id,tissue_color_hex)]
color=rbind(color,
	data.frame(Tissues=c('RPE (glu)','RPE (gal)'),
		tissue_site_detail_id=c('glucose','galactose'),
		tissue_color_hex=c('808080','D3D3D3')))
color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
mds=merge(mds,color[,list(tissue=tissue_site_detail_id,color=tissue_color_hex,`Tissues`)],by='tissue')
mds$Tissues=as.character(mds$Tissues)
mds[,label:=ifelse(`Tissues`%in%c('Liver','Muscle - Skeletal','Brain - Frontal Cortex (BA9)','Whole Blood','Cells - EBV-transformed lymphocytes','Pancreas','Testis','RPE (glu)'),`Tissues`,'')]
mds[,label:=str_replace(label,' - Frontal Cortex \\(BA9\\)','')]
mds[,label:=str_replace(label,'Cells - EBV-transformed lymphocytes','LCL')]
mds[,label:=str_replace(label,'RPE (glu)','RPE (glu)')]

# make scatter plot:
message('plotting...')
mds_centroid=unique(mds[,list(tissue,x_centroid,y_centroid,`Tissues`,label)])
color_map=color$tissue_color_hex
names(color_map)=color$`Tissues`
p1=ggplot(mds,aes(x=x,y=y,color=`Tissues`))+geom_point()+scale_color_manual(values=color_map)+geom_text(data=mds_centroid,aes(x=x_centroid,y=y_centroid,label=label),color='black')+theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+theme(axis.title=element_text(size=15))


mds[,label:=ifelse(Tissues%in%'RPE (glu)',sample,'')]
p2=ggplot(mds,aes(x=x,y=y,color=`Tissues`))+geom_point()+scale_color_manual(values=color_map)+geom_text_repel(aes(x=x,y=y,label=label),color='black')+theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+theme(axis.title=element_text(size=15))

pdf(paste0(figure_dir,'/mds.pdf'),height=6,width=14.5)
p1;p2
dev.off()


# Calculate PCs:
pc_res=prcomp(t(rpkm_mat))
saveRDS(list(pc_res,rpkm_mat,col_data,row_data),paste0(out_dir,'/pc.rds'))
pc=data.table(pc_res$x[,1:3],keep.rownames=T)

colnames(pc)=c('sample','PC1','PC2','PC3')
pc$tissue=col_data$tissue
# mds[,x_centroid:=median(x),by='tissue']
# mds[,y_centroid:=median(y),by='tissue']
pc[,PC1_centroid:=mean(PC1),by='tissue']
pc[,PC2_centroid:=mean(PC2),by='tissue']
pc[,PC3_centroid:=mean(PC3),by='tissue']


# prepare color palette: 
color=fread(color_file)[,list(`Tissues`=tissue_site_detail,tissue_site_detail_id,tissue_color_hex)]
color=rbind(color,
	data.frame(Tissues=c('RPE (glu)','RPE (gal)'),
		tissue_site_detail_id=c('glucose','galactose'),
		tissue_color_hex=c('808080','D3D3D3')))
color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
pc=merge(pc,color[,list(tissue=tissue_site_detail_id,color=tissue_color_hex,`Tissues`)],by='tissue')
pc$Tissues=as.character(pc$Tissues)
pc[,label:=ifelse(`Tissues`%in%c('Liver','Muscle - Skeletal','Brain - Frontal Cortex (BA9)','Whole Blood','Cells - EBV-transformed lymphocytes','Pancreas','Testis','RPE (glu)'),`Tissues`,'')]
pc[,label:=str_replace(label,' - Frontal Cortex \\(BA9\\)','')]
pc[,label:=str_replace(label,'Cells - EBV-transformed lymphocytes','LCL')]
pc[,label:=str_replace(label,'RPE (glu)','RPE (glu)')]

# make scatter plot:
message('plotting...')
pc_centroid=unique(pc[,list(tissue,x_centroid,y_centroid,`Tissues`,label)])
color_map=color$tissue_color_hex
names(color_map)=color$`Tissues`
p3=ggplot(pc,aes(x=PC1,y=PC2,color=`Tissues`))+geom_point()+scale_color_manual(values=color_map)+geom_text(data=mds_centroid,aes(x=x_centroid,y=y_centroid,label=label),color='black')+theme_bw()+xlab('PC1')+ylab('PC2')+theme(axis.title=element_text(size=15))
pc[,label:=ifelse(Tissues%in%'RPE (glu)',sample,'')]
p4=ggplot(pc,aes(x=PC1,y=PC2,color=`Tissues`,label=label))+geom_point()+scale_color_manual(values=color_map)+geom_text_repel(color='black')+theme_bw()+xlab('MDS Coordinate 1')+ylab('MDS Coordinate 2')+theme(axis.title=element_text(size=15))

pdf(paste0(figure_dir,'/pc.pdf'),width=14.5,height=6)
p3;p4
dev.off()
