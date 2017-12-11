library(data.table)
library(stringr)
library(ggplot2)

rpkm_fn='../data/rnaseq/rpkm/rpe.filt.rpkm'
rpe_signature_fn='rpe_marker_gene/rpe_signature.txt'
gencode_fn='/mnt/lab_data/montgomery/shared/annotations/gencode.v19.annotation.gtf'
fig_dir='../figures/rpe_marker_gene/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

select_treatment=function(rpkm,treatment,format=c('data.table','matrix')){
	
	if ('gene_id'%in%colnames(rpkm)){
		gene_id=rpkm[,list(gene_id)]
	}
	
	x=rpkm[,str_detect(colnames(rpkm),tolower(treatment)),with=FALSE]
	
	format=match.arg(format)
	if (format=='data.table'){
		cbind(gene_id,x)
	} else {
		x=as.matrix(x)
		rownames(x)=gene_id$gene_id
		return(x)
	}
	
}

get_gene_name_and_id=function(fn){
	x=fread(fn)
	x[,gene_id:=str_extract(V9,'(?<=gene_id ")(.+?)(?=";)')]
	x[,gene_name:=str_extract(V9,'(?<=gene_name ")(.+?)(?=";)')]
	y=unique(x[,list(gene_id,gene_name)])
	return(y)
}


read_rpe_signature=function(fn,gencode_fn){
	x=fread(fn)
	gencode=get_gene_name_and_id(gencode_fn)
	gene_id=gencode[gene_name%in%x$`Gene Symbol`,gene_id]
	return(gene_id)
}



plot_distn=function(data,title){
	p=ggplot(data,aes(x=rpkm,fill=type,color=type))+
	geom_histogram(aes(y=..density..),position=position_dodge(width=0.5),bins=20)+
	scale_x_log10()+theme_classic()+
	scale_fill_discrete(name='Gene',labels=c('Others','RPE marker'))+
	scale_color_discrete(guide='none')+
	xlab('RPKM')+ylab('Density')+ggtitle(title)
	return(p)
}




main=function(){
	rpkm=fread(rpkm_fn)
	rpe_signature=read_rpe_signature(rpe_signature_fn,gencode_fn)

	for(treatment in c('glucose','galactose')){
		x=select_treatment(rpkm,treatment,format='matrix')
		x_median=apply(x,1,median)
		data=data.table(rpkm=x_median,gene_id=names(x_median))
		data[,type:=ifelse(gene_id%in%rpe_signature,'RPE','Others')]

		p=plot_distn(data,treatment)
		pdf(sprintf('%s/rpkm_dist_%s.pdf',fig_dir,treatment),height=4,width=4)
		print(p)
		dev.off()
	}
}

main()