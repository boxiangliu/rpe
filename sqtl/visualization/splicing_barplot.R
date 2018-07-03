# Author: Mike Gloudemans
# 6/21/2018
# Plot splicing differences among individuals at RDH5 locus

library(cowplot)
library(data.table)
fig_dir = '../figures/sqtl/visualization/splicing_barplot/'
meta_fn='../data/meta/dna2rna.txt'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


key_events = c("chr12:56115278:56115473:clu_1202",
	       "chr12:56115278:56117670:clu_1202",
	       "chr12:56115731:56117670:clu_1202")
top_samps = c("X021512",
		"X041212.9319",
		"X070810",
		"X080410")

# Read data and subset down to the events we care about
d = read.table("/users/mgloud/projects/brain_gwas/scripts/auxiliary/rpe/splicing_plot/key_snps.txt", row.names=1, header=TRUE)
d = d[rownames(d) %in% key_events,]

# Compute splicing proportions relative to total number of events
d = t(t(d) / colSums(d))

# Sort samples first by genotype, then alphabetically
d = d[, order(colnames(d))]
d = d[, order(colnames(d) %in% top_samps)]
colnames(d) = gsub("\\.", "-", colnames(d))
colnames(d) = gsub("X", "", colnames(d))

# Assemble into data frame for plotting
new_data = data.frame(skipping = d[2,], sample=colnames(d))
new_data$sample = factor(new_data$sample, levels=new_data$sample)
new_data$genotype = "CC"
new_data$genotype[20:23] = "CA"

# Update sample names:
meta = fread(meta_fn,header=TRUE,colClasses = c('character','character','integer'))
new_data = merge(new_data,meta[,list(RNA,NAME)],by.x='sample',by.y='RNA')

setorder(new_data,-genotype,NAME)
new_data$NAME = factor(new_data$NAME,levels=new_data$NAME)

# Plot
p = ggplot(new_data, aes(NAME, skipping, fill=genotype)) +
	geom_col(color="black") +
	theme(axis.text.x=element_text(angle=45, hjust=1)) +
	scale_x_discrete(name="Sample ID") +
	scale_y_continuous(name="Proportion exon skipping") + 
	theme(legend.position = c(0.05, 0.8)) +
	xlab(NULL) + 
	scale_fill_discrete(name = 'rs3138141')

fig_fn = sprintf('%s/splicing_barplot.pdf',fig_dir)
ggsave(fig_fn, width=8, height=4)

