# Normalize using the GTEx pipeline. 
library(preprocessCore)

# Functions: 
getInverseNormal <- function(x){
        x <- as.numeric(x)
        xx2<-qnorm((rank(x,na.last = "keep") - 0.5) / sum(!is.na(x)))
        return(xx2)
}


# Command args: 
in_f='../data/rnaseq/rpkm/rpe.filt.rpkm'
out_f='../data/rnaseq/rpkm/rpe.filt.quant_norm.rpkm'


# Main:
rpkm=read.table(in_f,header=T,check.names=F)
rownames(rpkm)=rpkm[,1]
rpkm[,1]=NULL
rpkm_qnorm=as.data.frame(normalize.quantiles(as.matrix(rpkm)))
colnames(rpkm_qnorm)=colnames(rpkm)
rpkm_qnorm$gene_id=rownames(rpkm)
rpkm_qnorm=rpkm_qnorm[,c(ncol(rpkm_qnorm),1:(ncol(rpkm_qnorm)-1))]
write.table(rpkm_qnorm,out_f,quote=F,sep='\t',row.names=F,col.names=T)