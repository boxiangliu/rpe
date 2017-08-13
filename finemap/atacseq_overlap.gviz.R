library(data.table)
library(stringr)
library(cowplot)
library(rtracklayer)
library(ggrepel)
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)


bigwig_fn=c("../data/atacseq/rpe/out/signal/macs2/rep1/Sample1_RPE_R1.trim.PE2SE.nodup.tn5.pf.fc.signal.bigwig")
dbsnp_fn='../data/dbsnp/snp147.txt'
fig_dir='../figures/finemap/atacseq_overlap/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

# Function: 
get_bw=function(f,seq,start,end){
	gr=GRanges(seqnames = seq,
		ranges = IRanges(start = start, end = end))

	bw=import(f,which=gr)
	bw=as.data.frame(bw)
	setDT(bw)
	
	temp1=bw[,list(seqnames,start,score)]
	temp2=bw[,list(seqnames,end,score)]
	setnames(temp1,'start','pos')
	setnames(temp2,'end','pos')
	bw2=rbind(temp1,temp2)

	setorder(bw2,seqnames,pos)
	return(bw2)
}


snp2GRanges=function(file){
	snp=fread(file,sep='\t',skip=1,select=c(2,3,4,5),
		col.names=c('chr','start','end','name'))
	snp=GenomicRanges::makeGRangesFromDataFrame(snp,keep.extra.columns=TRUE)
	return(snp)
}


combine_data=function(bigwig_fn,gwas,dbsnp_fn){
	gen='hg19'

	tracks=list()

	bw_track = DataTrack(range = bigwig_fn, genome = gen, 
		window = -1, name = 'ATAC', type='histogram',frame=TRUE,col.frame='black')

	tracks[['ATAC']]=bw_track


	print('gwas')
	gwas_track = DataTrack(data = gwas$logp, chromosome = gwas$chr, start = gwas$pos, 
		end = gwas$pos, id = gwas$rsid, genome = gen, name = 'GWAS', type='p',frame=TRUE,col.frame='black')
	tracks[['gwas']] = gwas_track


	# print('chromHMM')
	# chromHMM=chromHMM2GRanges(chromHMM_fn)
	# ch_track = AnnotationTrack(range = chromHMM, genome = gen, 
	# 	name = 'ChromHMM', feature = chromHMM$itemRgb, id = chromHMM$name, 
	# 	stacking = 'dense', collapse = FALSE, frame=TRUE,col.frame='black',groupAnnotation = 'id')

	# feat = unique(feature(ch_track))
	# featCol = setNames(as.list(rgb(t(sapply(strsplit(feat, ","),
	# 	as.numeric)), maxColorValue=255)), feat)
	# displayPars(ch_track) = featCol

	# tracks[['chromHMM']] = ch_track

	print('Gene model')
	grtrack = GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, genome = gen,name = 'Gene Model',  transcriptAnnotation = "symbol",frame=TRUE,col.frame='black')
	symbols = unlist(mapIds(org.Hs.eg.db, gene(grtrack), "SYMBOL", "ENTREZID", multiVals = "first"))
	symbol(grtrack) = symbols[gene(grtrack)]
	tracks[['gene_model']] = grtrack


	print('dbSNP')
	snp=snp2GRanges(dbsnp_fn)
	snp_track = AnnotationTrack(range = snp, genome = gen, name = 'SNPs', id = snp$name, 
		stacking = 'full', collapse = FALSE, frame=TRUE,col.frame='black', groupAnnotation = 'id')
	tracks[['snp']] = snp_track


	print('Genomic location')
	gtrack = GenomeAxisTrack(labelPos = "below")
	tracks[['genomic_location']] = gtrack

	return(tracks)
}


plot_ase=function(x,title){
	pi=x[,pi]
	freq=x[,freq]
	ref=x[,ref]
	alt=x[,alt]

	to_plot=data.frame(allele=c(sprintf('ref: %s (%.02f)',ref,1-freq),sprintf('alt: %s (%.02f)',alt,freq)),ase=c(1-pi,pi))
	p=ggplot(to_plot,aes(allele,ase,label=sprintf('%.02f%%',ase*100)))+geom_bar(stat='identity')+geom_text(nudge_y=0.02)+ggtitle(title)
	return(p)
}


# Read GWAS data:
gwas_fn='../data/gwas/Fritsche_2015_AdvancedAMD.txt'
gwas=fread(gwas_fn,select=c(1:3,8),col.names=c('rsid','chr','pos','pval'))
gwas[,chr:=paste0('chr',chr)]
gwas[,logp:=-log10(pval)]


# Combine all data:
tracks=combine_data(bigwig_fn,gwas,dbsnp_fn)



# TCF21:
chr='chr12'
eqtl_fn='../processed_data/rasqual/output/glucose/joint/chr12/ENSG00000135437.5_RDH5.txt'


eqtl=fread(eqtl_fn,select=c(1:4,7,11,25),col.names=c('fid','sid','chr','pos','maf','chisq','r2_rsnp'))
eqtl[,pval:=pchisq(chisq,df=1,lower.tail=FALSE)]
eqtl[,logp:=-log10(pval)]

eqtl_track = DataTrack(data = eqtl$logp, start = eqtl$pos, 
	end = eqtl$pos, chromosome = chr, genome = 'hg19', name = 'RAS', type='p',frame=TRUE,col.frame='black')

tracks[['RASQUAL']]=eqtl_track
tracks=tracks[c(1,2,6,3:5)]

pdf(sprintf('%s/rdh5.pdf',fig_dir,chr,start,end),height=10,width=8)
start=5.56e7
end=5.66e7
plotTracks(tracks[-5], chromosome = chr, from = start, to = end)

# rs3138141 is 50bp downstream of the splice donor site.
start=5.611e7
end=5.612e7
ht = HighlightTrack(trackList=tracks, start = c(56115585,56115778)-25, width = 50, chromosome= chr)
plotTracks(list(ht), chromosome = chr, from = start, to = end,showFeatureId=FALSE,showId=TRUE)

start=5.621e7
end=5.622e7
ht = HighlightTrack(trackList=tracks, start = 56213297-25, width = 50, chromosome= chr)
plotTracks(list(ht), chromosome = chr, from = start, to = end,showFeatureId=FALSE,showId=TRUE)

start=5.6294e7
end=5.6352e7
ht = HighlightTrack(trackList=tracks, start = c(56294785,56329537,56335150,56351346)-25, width = 50, chromosome= chr)
plotTracks(list(ht), chromosome = chr, from = start, to = end,showFeatureId=FALSE,showId=TRUE)
dev.off()
