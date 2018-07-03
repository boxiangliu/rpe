library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)
library(manhattan)
library(cowplot)
library(ggrepel)
library(RCircos)
source('utils/genome_annotation.R')

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
treeQTL_MT_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
out_dir = '../processed_data/QTL_landscape/eQTL_manhattan/'
fig_dir = '../figures/QTL_landscape/eQTL_manhattan/'
imprinted_gene_fn = '../data/imprinted_genes/imprinted_genes.txt'

if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_rasqual = function(fn){
	rasqual = fread(fn,header=FALSE)[,c(1,3,4,11,12)]
	colnames = c('fid','chr','pos','chisq','pi')
	setnames(rasqual,colnames)
	rasqual[,rank := rank(-chisq,ties.method='random')]
	rasqual = rasqual[rank == 1]
	rasqual$rank = NULL
	split_fid = str_split_fixed(rasqual$fid,'_',2)
	rasqual$gene_id = split_fid[,1]
	rasqual$gene_name = split_fid[,2]
	rasqual$fid = NULL
	rasqual[,pval := pchisq(chisq,df=1,lower.tail=FALSE)]
	rasqual[,logp := -log10(pval)]
	return(rasqual)
}

read_TreeQTL_MT = function(fn){
	treeQTL_MT = fread(fn)
	return(treeQTL_MT)
}

add_eQTL_indicator = function(gene_id, eqtl_list){
	eqtl_indicator = gene_id %in% eqtl_list
	return(eqtl_indicator)
}

rbind_eQTLs = function(x,y){
	z = copy(y)
	z[,logp := -logp]
	rbind(x,z)
}

plot_manhattan = function(eqtl){
	data = copy(eqtl)
	build = 'hg19'
	data[,y:=logp]
	data[,chrom:=chr]
	data = add_cumulative_pos(data, build)
	color1 = 'black'
	color2 = 'grey'
	data = add_color(data, color1 = color1, color2 = color2)
	data[eQTL==TRUE,color := 'red']
	data[treatment_specific == TRUE,color:='blue']
	chrom_lengths = get_chrom_lengths(build)
	xmax = get_total_length(chrom_lengths)
	x_breaks = get_x_breaks(chrom_lengths)
	color_map = unique(data$color)
	names(color_map) = unique(data$color)

	p=ggplot(data, aes(x = cumulative_pos, y = y, color = color)) + 
		geom_point(alpha = 0.5) + 
		theme_classic() + 
		scale_x_continuous(limits = c(0, xmax), expand = c(0.01, 0), breaks = x_breaks, labels = names(x_breaks), name = "Chromosome") + 
		scale_y_continuous(expand = c(0.01, 0), name = expression("-log10(P-value)"),labels = abs) + 
		scale_color_manual(values = color_map, name = '', breaks = c('red','blue'), labels = c('eQTL','Treatment-specific eQTL')) + 
		theme(legend.position = 'top') +
		facet_grid(treatment~., scales = 'free_y') 
	return(p)
}

label_eQTL_by_logp = function(data,top=10){
	diff = data[glucose_specific==TRUE,list(gene_id, gene_name, diff = glucose_logp - galactose_logp)]
	setorder(diff,-diff)
	data$label = ''
	data[,label := ifelse(gene_name %in% diff[1:top,gene_name], gene_name, label)]
	
	diff = data[galactose_specific==TRUE,list(gene_id, gene_name, diff = galactose_logp - glucose_logp)]
	setorder(diff,-diff)
	data[,label := ifelse(gene_name %in% diff[1:top,gene_name], gene_name, label)]

	sum = data[galactose_shared&glucose_shared,list(gene_id, gene_name, sum = galactose_logp + glucose_logp)]
	setorder(sum,-sum)
	data[,label := ifelse(gene_name %in% sum[1:top,gene_name], gene_name, label)]
	return(data)
}

merge_glucose_galactose = function(glucose,galactose){
	data = merge(
		glucose[,list(gene_id,gene_name,glucose_pi = pi, glucose_logp=logp,glucose_specific=treatment_specific,glucose_shared = eQTL&!treatment_specific)],
		galactose[,list(gene_id,galactose_pi = pi, galactose_logp=logp,galactose_specific=treatment_specific,galactose_shared = eQTL&!treatment_specific)]
		,by='gene_id')
	return(data)
}

rank_eQTL_by_pi = function(data){
	delta_pi = abs(data$glucose_pi - 0.5) - abs(data$galactose_pi - 0.5)
	delta_pi = ifelse(data$glucose_specific==TRUE,delta_pi,-Inf)
	data$glucose_specific_rank = rank(-delta_pi,ties.method='first')

	delta_pi = abs(data$galactose_pi - 0.5) - abs(data$glucose_pi - 0.5)
	delta_pi = ifelse(data$galactose_specific==TRUE,delta_pi,-Inf)
	data$galactose_specific_rank = rank(-delta_pi,ties.method='first')

	sum_pi = abs(data$glucose_pi - 0.5) + abs(data$galactose_pi - 0.5)
	sum_pi = ifelse(data$glucose_shared==TRUE,sum_pi, -Inf)
	data$shared_rank = rank(-sum_pi,ties.method='first')
	return(data)
}

plot_scatterplot = function(glucose,galactose,top = 5, log = TRUE){
	data = merge_glucose_galactose(glucose,galactose)
	imprinted_gene = fread(imprinted_gene_fn)
	data = data[!(gene_name %in% imprinted_gene$Gene)]
	data = label_eQTL_by_logp(data,top = top)
	data$color = 'black'
	data[,color := ifelse(glucose_shared,'green',color)]
	data[,color := ifelse(galactose_shared,'green',color)]
	data[,color := ifelse(glucose_specific,'red',color)]
	data[,color := ifelse(galactose_specific,'blue',color)]
	if (log){
		p = ggplot(data,aes(glucose_logp,galactose_logp,color=color,label=label)) + 
			geom_point(alpha = 0.5) + 
			geom_abline(intercept = 0, slope =1, color = 'black', linetype = 'dashed') + 
			scale_color_manual(name = '', breaks = c('blue','red','green'), labels = c('Galactose-specific','Glucose-specific','Shared'),values = c(black='black',red='red',blue='blue',green='green')) + 
			xlab(expression(Glucose -log['10'](P-value))) + 
			ylab(expression(Galactose -log['10'](P-value))) + 
			theme(legend.position = 'top') + 
			geom_text_repel(color = 'black', show_guide =FALSE) + 
			scale_x_log10() + 
			scale_y_log10()
	} else {
		p = ggplot(data,aes(glucose_logp,galactose_logp,color=color,label=label)) + 
			geom_point(alpha = 0.5) + 
			geom_abline(intercept = 0, slope =1, color = 'black', linetype = 'dashed') + 
			scale_color_manual(name = '', breaks = c('blue','red','green'), labels = c('Galactose-specific','Glucose-specific','Shared'),values = c(black='black',red='red',blue='blue',green='green')) + 
			xlab(expression(Glucose -log['10'](P-value))) + 
			ylab(expression(Galactose -log['10'](P-value))) + 
			theme(legend.position = 'top') + 
			geom_text_repel(color = 'black',show_guide = FALSE)
	}


	return(p)
}

plot_scatterplot_label_by_pi = function(glucose,galactose,top = 5, log = TRUE){
	data = merge_glucose_galactose(glucose,galactose)
	imprinted_gene = fread(imprinted_gene_fn)
	data = data[!(gene_name %in% imprinted_gene$Gene)]
	data = rank_eQTL_by_pi(data)
	data[,label := ifelse(glucose_specific_rank<=5|galactose_specific_rank<=5|shared_rank<=5,gene_name,'')]
	data$color = 'black'
	data[,color := ifelse(glucose_shared,'green',color)]
	data[,color := ifelse(galactose_shared,'green',color)]
	data[,color := ifelse(glucose_specific,'red',color)]
	data[,color := ifelse(galactose_specific,'blue',color)]
	if (log){
		p = ggplot(data,aes(glucose_logp,galactose_logp,color=color,label=label)) + 
			geom_point(alpha = 0.5) + 
			geom_abline(intercept = 0, slope =1, color = 'black', linetype = 'dashed') + 
			scale_color_manual(name = '', breaks = c('blue','red','green'), labels = c('Galactose-specific','Glucose-specific','Shared'),values = c(black='black',red='red',blue='blue',green='green')) + 
			xlab(expression(Glucose -log['10'](P-value))) + 
			ylab(expression(Galactose -log['10'](P-value))) + 
			theme(legend.position = 'top') + 
			geom_text_repel(color = 'black', show.legend =FALSE) + 
			scale_x_log10() + 
			scale_y_log10()
	} else {
		p = ggplot(data,aes(glucose_logp,galactose_logp,color=color,label=label)) + 
			geom_point(alpha = 0.5) + 
			geom_abline(intercept = 0, slope =1, color = 'black', linetype = 'dashed') + 
			scale_color_manual(name = '', breaks = c('blue','red','green'), labels = c('Galactose-specific','Glucose-specific','Shared'),values = c(black='black',red='red',blue='blue',green='green')) + 
			xlab(expression(Glucose -log['10'](P-value))) + 
			ylab(expression(Galactose -log['10'](P-value))) + 
			theme(legend.position = 'top') + 
			geom_text_repel(color = 'black',show.legend = FALSE)
	}
	return(p)
}

RCircos.Vertical.Line.Plot <- function(line.data=NULL, track.num=NULL, 
		side=c("in", "out"), line.width=1, inside.pos=NULL, outside.pos=NULL,
		genomic.columns=3, is.sorted=TRUE, color = 'black'){
	if(is.null(line.data)) 
		stop("Genomic data missing in RCircos.Vertical.Line.Plot.\n");
	if(is.null(genomic.columns) || genomic.columns < 2) 
		stop("Missing number of columns for genomic position.\n");

	boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
								outside.pos, FALSE);
	outerPos <- boundary[1];
	innerPos  <- boundary[2];

	if(length(line.width) == 1)
		line.width <- rep(line.width, nrow(line.data));

	line.data <- RCircos.Get.Single.Point.Positions(line.data, 
						genomic.columns);
	point.index <- as.numeric(line.data$Location);
	
	RCircos.Pos <- RCircos.Get.Plot.Positions();
	RCircos.Par <- RCircos.Get.Plot.Parameters();
	line.colors <- RCircos.Get.Plot.Colors(line.data, color);

	for(a.point in seq_len((nrow(line.data))))
	{
		lines(  c(RCircos.Pos[point.index[a.point] ,1]*innerPos,
					RCircos.Pos[point.index[a.point], 1]*outerPos),
				c(  RCircos.Pos[point.index[a.point], 2]*innerPos,
					RCircos.Pos[point.index[a.point], 2]*outerPos),
				col=line.colors[a.point], lwd=line.width[a.point]);
	}
}


setup_circos = function(tracks.inside, tracks.outside, chr.exclude = NULL){
	data(UCSC.HG19.Human.CytoBandIdeogram)
	cyto.info = UCSC.HG19.Human.CytoBandIdeogram
	RCircos.Set.Core.Components(cyto.info = cyto.info, chr.exclude = chr.exclude, tracks.inside = tracks.inside, tracks.outside = tracks.outside)
}

RCircos.Label.Chromosome.Names = function(chr.name.pos = NULL){
	RCircos.Par <- RCircos.Get.Plot.Parameters()
	if (is.null(chr.name.pos)) {
		chr.name.pos <- RCircos.Par$chr.name.pos
	}
	else {
		if (chr.name.pos < RCircos.Par$track.in.start) 
			stop("Chromosome name positions overlap with inside track.\n")
		if (chr.name.pos > RCircos.Par$track.out.start) 
			message("May not plot data tracks outside of chromosomes.\n")
	}
	RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
	RCircos.Pos <- RCircos.Get.Plot.Positions()
	chroms <- unique(RCircos.Cyto$Chromosome)
	rightSide <- nrow(RCircos.Pos)/2
	for (aChr in seq_len(length(chroms))) {
		chrRows <- which(RCircos.Cyto$Chromosome == chroms[aChr])
		chrStart <- min(RCircos.Cyto$StartPoint[chrRows])
		chrEnd <- max(RCircos.Cyto$EndPoint[chrRows])
		mid <- round((chrEnd - chrStart + 1)/2, digits = 0) + 
			chrStart
		chrName <- sub(pattern = "chr", replacement = "", chroms[aChr])
		text(RCircos.Pos[mid, 1] * chr.name.pos, 
			RCircos.Pos[mid, 2] * chr.name.pos, label = chrName, 
			pos = ifelse(mid <= rightSide, 4, 2), srt = RCircos.Pos$degree[mid])
	}
}

plot_circos = function(glucose,galactose,rpe_specific_eGenes){

	RCircos.Set.Plot.Area()
	RCircos.Highligh.Chromosome.Ideogram(highlight.pos = 1.47)
	RCircos.Label.Chromosome.Names(chr.name.pos = 1.49)
	# RCircos.Chromosome.Ideogram.Plot()

	data = merge(
		glucose[,list(gene_id, chr, start = pos, end = pos, gene_name, glucose = eQTL)],
		galactose[,list(gene_id,galactose = eQTL)],
		by='gene_id')
	data[,shared := glucose & galactose]
	glucose_line = data[glucose == TRUE & shared == FALSE, list(chr,start,end,glucose)]
	galactose_line = data[galactose == TRUE & shared == FALSE, list(chr,start,end,galactose)]
	shared_line = data[shared == TRUE, list(chr,start,end,shared)]
	rpe_specific_line = data[gene_id %in% rpe_specific_eGenes$gene_id,list(chr,start,end,gene_name)]
	RCircos.Vertical.Line.Plot(shared_line, track.num=1, side="in",color='black')
	RCircos.Vertical.Line.Plot(glucose_line, track.num=2, side="in",color='#F87660')
	RCircos.Vertical.Line.Plot(galactose_line, track.num=3, side="in",color='#619CFF')
	RCircos.Vertical.Line.Plot(rpe_specific_line, track.num=4, side="in",color='#C77CFF')
	setDF(rpe_specific_line)
	RCircos.Gene.Name.Plot(gene.data = rpe_specific_line,name.col=4,track.num=5, side='in')
}

get_response_eQTLs = function(glucose,galactose){
	data = merge(
		glucose[,list(gene_id,gene_name,glucose_logp=logp,glucose_specific=treatment_specific,glucose_shared = eQTL)],
		galactose[,list(gene_id,galactose_logp=logp,galactose_specific=treatment_specific,galactose_shared = eQTL)]
		,by='gene_id')
	glucose_specific = data[glucose_specific==TRUE,list(gene_id, gene_name, metric = glucose_logp - galactose_logp)]
	setorder(glucose_specific,-metric)
	glucose_specific$eQTL = 'Glucose-specific'
	
	galactose_specific = data[galactose_specific==TRUE,list(gene_id, gene_name, metric = galactose_logp - glucose_logp)]
	setorder(galactose_specific,-metric)
	galactose_specific$eQTL = 'Galactose-specific'

	shared = data[glucose_shared==TRUE&galactose_shared==TRUE,list(gene_id, gene_name, metric = galactose_logp + glucose_logp)]
	setorder(shared,-metric)
	shared$eQTL = 'Shared'
	response_eQTLs = rbind(glucose_specific,galactose_specific,shared)
	return(response_eQTLs)
}
# Main:
#-----------#
# read data #
#-----------#
treeQTL_MT = read_TreeQTL_MT(treeQTL_MT_fn)

glucose_list = list.files(glucose_dir,pattern = 'txt',recursive=TRUE,full.names=TRUE)
glucose = foreach(fn = glucose_list, .combine = rbind)%dopar%{
	rasqual = read_rasqual(fn)
	return(rasqual)
}

glucose$eQTL = add_eQTL_indicator(glucose$gene_id,treeQTL_MT[glucose==1,gene])
glucose$treatment_specific = add_eQTL_indicator(glucose$gene_id,treeQTL_MT[glucose==1&galactose==0,gene])
glucose$treatment = 'Glucose'
out_fn = sprintf('%s/rasqual_glucose_top_eQTL.txt',out_dir)
fwrite(glucose,out_fn,sep='\t')
glucose = fread(out_fn)

galactose_list = list.files(galactose_dir,pattern = 'txt',recursive=TRUE,full.names=TRUE)
galactose = foreach(fn = galactose_list, .combine = rbind)%dopar%{
	rasqual = read_rasqual(fn)
	return(rasqual)
}
galactose$eQTL = add_eQTL_indicator(galactose$gene_id,treeQTL_MT[galactose==1,gene])
galactose$treatment_specific = add_eQTL_indicator(galactose$gene_id,treeQTL_MT[galactose==1&glucose==0,gene])
galactose$treatment = 'Galactose'
out_fn = sprintf('%s/rasqual_galactose_top_eQTL.txt',out_dir)
fwrite(galactose,out_fn,sep='\t')
galactose = fread(out_fn)

#-----------#
# manhattan #
#-----------#
eqtl = rbind_eQTLs(galactose,glucose)
p = plot_manhattan(eqtl)
fig_fn = sprintf('%s/manhattan.pdf',fig_dir)
pdf(fig_fn, width = 8, height = 4)
print(p)
dev.off()

#--------------#
# scatter plot #
#--------------#
###--- ranked by difference in P-value ---###
p = plot_scatterplot(glucose,galactose,log=FALSE)
fig_fn = sprintf('%s/scatterplot.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)

p = plot_scatterplot(glucose,galactose,log=TRUE)
fig_fn = sprintf('%s/scatterplot_log.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)


###---- ranked by difference in pi ---####
p = plot_scatterplot_label_by_pi(glucose,galactose,log=FALSE)
fig_fn = sprintf('%s/scatterplot_label_by_pi.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)

p = plot_scatterplot_label_by_pi(glucose,galactose,log=TRUE)
fig_fn = sprintf('%s/scatterplot_log_label_by_pi.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)


#-------------#
# cicros plot #
#-------------#
rpe_specific_eGenes_fn = '../processed_data/rpe_specific_eQTL/specific_eGenes_v2//rpe_specific_eGenes.txt'
rpe_specific_eGenes = fread(rpe_specific_eGenes_fn)

mean_expression = read_mean_expression()
gene_annotation = read_gencode()
mean_expression = merge(mean_expression,gene_annotation,by=c('gene_id','gene_name'))
expressed_genes = mean_expression[mean_rpkm>0.5,]

num_shared_eQTL = nrow(glucose[eQTL == TRUE & treatment_specific==FALSE & gene_id %in% expressed_genes$gene_id])
num_glucose_eQTL = nrow(glucose[treatment_specific==TRUE & gene_id %in% expressed_genes$gene_id])
num_galactose_eQTL = nrow(galactose[treatment_specific==TRUE & gene_id %in% expressed_genes$gene_id])
num_rpe_specific_eQTL = length(unique(rpe_specific_eGenes$gene_id))
legend = c(paste0('Shared = ',num_shared_eQTL),
	paste0('Glucose = ',num_glucose_eQTL),
	paste0('Galactose = ',num_galactose_eQTL),
	paste0('RPE-specific = ',num_rpe_specific_eQTL))

setup_circos(tracks.inside = 6, tracks.outside = 0, chr.exclude = c('chrX','chrY'))
fig_fn = sprintf('%s/circos.pdf',fig_dir)
pdf(fig_fn,width = 6,height=6)
plot_circos(glucose,galactose,rpe_specific_eGenes)
legend(x = -0.7,y=0.5,bty = 'n', legend = legend,col = c('black','#F87660','#619CFF','#C77CFF'),pch = 19)
dev.off()