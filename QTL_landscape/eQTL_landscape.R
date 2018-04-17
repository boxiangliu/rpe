library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)
library(manhattan)
library(cowplot)
library(ggrepel)
library(RCircos)

glucose_dir = '../processed_data/rasqual/output/glucose/joint/'
galactose_dir = '../processed_data/rasqual/output/galactose/joint/'
treeQTL_MT_fn = '../processed_data/response_eQTL/treeQTL_MT/eGenesMT.txt'
out_dir = '../processed_data/QTL_landscape/eQTL_manhattan/'
fig_dir = '../figures/QTL_landscape/eQTL_manhattan/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

read_rasqual = function(fn){
	rasqual = fread(fn,header=FALSE)[,c(1,3,4,11)]
	colnames = c('fid','chr','pos','chisq')
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

label_eQTL = function(data,top=10){
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

plot_scatterplot = function(glucose,galactose,top = 5, log = TRUE){
	data = merge(
		glucose[,list(gene_id,gene_name,glucose_logp=logp,glucose_specific=treatment_specific,glucose_shared = eQTL)],
		galactose[,list(gene_id,galactose_logp=logp,galactose_specific=treatment_specific,galactose_shared = eQTL)]
		,by='gene_id')
	data = label_eQTL(data,top = top)
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
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
	RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)
}

plot_circos = function(glucose,galactose){
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot()
	data = merge(
		glucose[,list(gene_id, chr, start = pos, end = pos, gene_name, glucose = eQTL)],
		galactose[,list(gene_id,galactose = eQTL)],
		by='gene_id')
	data[,shared := glucose & galactose]
	glucose_line = data[glucose == TRUE, list(chr,start,end,glucose)]
	galactose_line = data[galactose == TRUE, list(chr,start,end,galactose)]
	shared_line = data[shared == TRUE, list(chr,start,end,shared)]
	RCircos.Vertical.Line.Plot(glucose_line, track.num=1, side="in",color='blue')
	RCircos.Vertical.Line.Plot(galactose_line, track.num=2, side="in",color='red')
	RCircos.Vertical.Line.Plot(shared_line, track.num=3, side="in",color='black')
}

get_response_eQTLs = function(glucose,galactose){
	data = merge(
		glucose[,list(gene_id,gene_name,glucose_logp=logp,glucose_specific=treatment_specific,glucose_shared = eQTL)],
		galactose[,list(gene_id,galactose_logp=logp,galactose_specific=treatment_specific,galactose_shared = eQTL)]
		,by='gene_id')
	glucose_specific = data[glucose_specific==TRUE,list(gene_id, gene_name, diff = glucose_logp - galactose_logp)]
	setorder(glucose_specific,-diff)
	glucose_specific$eQTL = 'Glucose-specific'
	
	galactose_specific = data[galactose_specific==TRUE,list(gene_id, gene_name, diff = galactose_logp - glucose_logp)]
	setorder(galactose_specific,-diff)
	galactose_specific$eQTL = 'Galactose-specific'

	response_eQTLs = rbind(glucose_specific,galactose_specific)
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

galactose_list = list.files(galactose_dir,pattern = 'txt',recursive=TRUE,full.names=TRUE)
galactose = foreach(fn = galactose_list, .combine = rbind)%dopar%{
	rasqual = read_rasqual(fn)
	return(rasqual)
}
galactose$eQTL = add_eQTL_indicator(galactose$gene_id,treeQTL_MT[galactose==1,gene])
galactose$treatment_specific = add_eQTL_indicator(galactose$gene_id,treeQTL_MT[galactose==1&glucose==0,gene])
galactose$treatment = 'Galactose'

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
p = plot_scatterplot(glucose,galactose,log=FALSE)
fig_fn = sprintf('%s/scatterplot.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)

p = plot_scatterplot(glucose,galactose,log=TRUE)
fig_fn = sprintf('%s/scatterplot_log.pdf',fig_dir)
save_plot(fig_fn,p,base_width = 6, base_height = 6)

# Save a list of response eQTLs:
response_eQTLs = get_response_eQTLs(glucose,galactose)
out_fn = sprintf('%s/response_eQTLs.txt',out_dir)
fwrite(response_eQTLs,out_fn,sep='\t')

#-------------#
# cicros plot #
#-------------#
setup_circos(tracks.inside = 3, tracks.outside = 0, chr.exclude = c('chrX','chrY'))
fig_fn = sprintf('%s/circos.pdf',fig_dir)
pdf(fig_fn,width = 6,height=6)
plot_circos(glucose,galactose)
legend('topright',legend = c('Glucose','Galactose','Shared'),col = c('blue','red','black'),pch = 19)
dev.off()