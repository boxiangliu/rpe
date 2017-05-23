#########
# setup: 
#########
library(reshape2)
library(stringr)
library(dplyr)
library(data.table)


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
args = commandArgs(TRUE)
message(sprintf('[ cmd ]: %s', args))
count_dir = sprintf('%s/',args[1])
figure_path = args[2]
print(count_dir)
print(figure_path)
filenames = list.files(count_dir, pattern = '*.Aligned.out.htseqcounts.txt', full.name = TRUE)

count = data.frame(
	gene_name = character(),
	count = integer(),
	sample_name = character(),
	stringsAsFactors = TRUE)

for (filename in filenames){
	temp = fread(filename) %>% 
	dplyr::rename(gene_name = V1, count = V2) %>% 
	mutate(sample_name = str_split_fixed(basename(filename),'\\.', n = 5)[1])
	count = rbind(count, temp)
}

count = dcast(count, gene_name ~ sample_name, value.var = 'count')
count = as.data.frame(count)

stopifnot(class(count) == 'data.frame')
count = count[rowSums(count[,2:ncol(count)]) != 0, ] # remove unexpressed genes. 
count = count[!(count$gene_name %in% c('__no_feature','__ambiguous','__alignment_not_unique',' __too_low_aQual','__not_aligned')),]
rownames(count) = count$gene_name
count$gene_name = NULL

# normalize (dividing by the column sum):
norm_count = scale(count, center = F, scale = colSums(count))

# quantile normalize: 
library(preprocessCore)
qnorm_count = normalize.quantiles(as.matrix(count))
colnames(qnorm_count) = colnames(count)

# log transform: 
log2_count = log2(count + 1/2) 

# regularized log transform: 
# library(DESeq2)
# sampleNames = str_extract(filenames, '[0-9]+_(Galactose|Glucose)')
# sampleCondition = str_extract(filenames, '(Galactose|Glucose)')
# sampleTable = data.frame(sampleName = sampleNames, fileName = filenames, condition = sampleCondition)
# dds = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = '.', design= ~ condition)
# dds = dds[ rowSums(counts(dds)) > 1, ] # filter out rows with total counts lower than 1
# rld = rlog(dds)

##save
# write.table(count, 'counts/htseqcounts.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

################################
# hierarchical clustering
# using spearman correlation and average linkage: 
################################
c = cor(qnorm_count, method="pearson")
d = as.dist(1-c)
hr = hclust(d, method = "average", members=NULL) 
pdf(figure_path)
par(mar=c(7.1,4.1,4.1,2.1))
plot(as.dendrogram(hr), edgePar=list(col=3, lwd=4), horiz=FALSE, main = 'hierarchical clustering: pearson correlation, average linkage') 

##############################################
# pairwise scatterplot of log10(gene counts) of 5 samples
# put pearson correlation on upper panel
##############################################
par(mar = c(5, 4, 4, 2) + 0.1)
pairs(log2_count, upper.panel = panel.cor, main = 'log2(gene_count + 1/2)')

#######################
# rank-frequency plot
#######################
sorted_count = count
for (i in 1:ncol(count)){sorted_count[,i] = sort(count[,i],decreasing=T)}
plot(1:nrow(sorted_count), sorted_count[,1], col = 1, log='y', ylim = c(1,460000), xlab = 'rank', ylab = 'count')
points(1:nrow(sorted_count), sorted_count[,2], col = 2)
points(1:nrow(sorted_count), sorted_count[,3], col = 3)
points(1:nrow(sorted_count), sorted_count[,4], col = 4)
legend('topright',legend = colnames(sorted_count), pch = 1, col = c(1,2,3,4))
dev.off()
