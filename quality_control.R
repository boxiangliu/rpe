#########
# setup: 
#########
setwd('/srv/gs1/projects/montgomery/bliu2/rpe_rna_seq')
library(reshape2)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gplots)

#############
# FUNCTIONS
#############

## pairs function add-on. 
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
 usr <- par("usr"); on.exit(par(usr))
 par(usr = c(0, 1, 0, 1))
 r <- abs(cor(x, y))
 txt <- format(c(r, 0.123456789), digits = digits)[1]
 txt <- paste0(prefix, txt)
 if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
 text(0.5, 0.5, txt, cex = cex.cor * r)
}

######################
# read in count data:
######################
filenames = list.files('data/count/', pattern = '*.merged.Aligned.out.sorted.count', full.name = TRUE)
count = data.frame(
	gene_name = character(),
	count = integer(),
	sample_name = character(),
	stringsAsFactors = TRUE)

for (filename in filenames){
	temp = fread(filename) %>% 
	dplyr::rename(gene_name = V1, count = V2) %>% 
	mutate(sample_name = str_split_fixed(basename(filename),'_', n = 2)[1])
	count = rbind(count, temp)
}

count = dcast(count, gene_name ~ sample_name, value.var = 'count')
count = count[rowSums(count[2:ncol(count)]) != 0, ] # remove unexpressed genes. 

##save
write.table(count, 'processed_data/count_rna_seq_rpe.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

################################
# hierarchical clustering
# using spearman correlation and average linkage: 
################################

c = cor(count[2:6], method="spearman")
d = as.dist(1-c)
hr = hclust(d, method = "average", members=NULL) 
pdf('figures/quality_control.R-hierarchical_clustering-150318.pdf')
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=FALSE, main = 'hierarchical clustering: pearson correlation, average linkage') 
dev.off()

#################################################################
# scatterplot dARPE19 with ARPE_DMSO gene counts (control sample): 
#################################################################

dARPE19 = fread('data/count/dARPE19_S5.merged.Aligned.out.sorted.count') %>% 
dplyr::rename(gene_name = V1, count = V2) %>% 
slice(1:(nrow(.)-5)) %>% 
filter(count != 0)

ARPE_DMSO = fread('data/ARPE_control_sample/ARPE-19-DMSO-2.txt')
to_plot = merge(dARPE19, ARPE_DMSO, by = 'gene_name') %>%
select(
	gene_name,
	dARPE19_count = count,
	ARPE_DMSO_count = read_count)

pdf('figures/quality_control.R-gene_count_dARPE19_vs_ARPE_DMSO-150315.pdf')
qplot(x = ARPE_DMSO_count, y = dARPE19_count, data = to_plot, log = 'xy') + 
geom_abline(intercept = 0, slope = 1, color = 'red') +
ggtitle('gene count: dARPE19 and ARPE_DMSO') + 
xlab('ARPE DMSO') +
ylab('dARPE19') 
dev.off()

#################################################
# plot gene expression level in descending order: 
#################################################

dARPE19 = dARPE19  %>% arrange(desc(count))
ARPE_DMSO = ARPE_DMSO %>% arrange(desc(read_count))

pdf('figures/quality_control.R-gene_count_vs_index_dARPE19_vs_ARPE_DMSO-150315.pdf')
qplot(x = 1:nrow(ARPE_DMSO), y = ARPE_DMSO$read_count, log = 'y', color = 'blue') + 
geom_point(aes(x = 1:nrow(dARPE19), y = dARPE19$count, color = 'red')) +
xlab('index') + 
ylab('gene count') + 
scale_color_discrete(name="sample",
labels=c("ARPE_DMSO", "dARPE19"))
dev.off()

##############################################
# pairwise scatterplot of log10(gene counts) of 5 samples
# put pearson correlation on upper panel
##############################################
log10_count = log10(count[,2:6] + 0.1)
pdf('figures/quality_control.R-pairwise_scatterplot_gene_count_all_samples.pdf')
pairs(log10_count, upper.panel = panel.cor)
dev.off()
