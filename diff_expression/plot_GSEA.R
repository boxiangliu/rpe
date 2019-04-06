library(cowplot)
library(Hmisc)
library(stringr)
library(data.table)

fig_dir = '../figures/diff_expression/plot_GSEA/'
out_dir = '../processed_data/diff_expression/plot_GSEA/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

# Read in data
data = read.table("diff_expression/RPEgsea_summary.txt", header = T)

# Set FDR.q.vals of zero to 1 divided by the number of permutations
data[838,6] = 1/10000
data[846,6] = 1/10000


beautify_pathway_names = function(name){
	name = gsub('_',' ',name)
	split_name = str_split_fixed(name,' ',2)
	split_name[,1] = capitalize(tolower(split_name[,1]))
	split_name[,2] = capitalize(tolower(split_name[,2]))
	name = split_name[,2]
	split_name = str_split_fixed(name,' ',2)
	name = paste(split_name[,1],split_name[,2],sep='\n')
	name = str_replace(name,'[rR]rna','rRNA')
	name = str_replace(name,'[Mm]torc1','mTORC1')
	name = str_replace(name,'Hypoxia\n','Hypoxia')
	name = str_replace(name,'rRNA\nmetabolic process','rRNA metabolic\nprocess')
	# name = str_replace(name,'\\(Go\\)','\\(GO\\)')
	# name = str_replace(name,'Ribonucleoprotein complex biogenesis','Ribonucleoprotein\ncomplex biogenesis')
	return(name)
}

data$NAME = beautify_pathway_names(data$NAME)


# Recover significant pathways and plot them
data_sig = data[which(data$FDR.q.val < 0.05),]
data_plot = data_sig[,c(1,2,4,6,12)]
setDT(data_plot)
names(data_plot) = c("Pathway","Size","NES","FDR","Direction")
setorder(data_plot,-capitalize(Pathway))
fwrite(data_plot,'../processed_data/diff_expression/plot_GSEA/data_plot.txt', sep = '\t')

data_plot = data_plot[rev(order(capitalize(data_plot$Pathway))),]
data_plot[,Pathway:=factor(Pathway,levels=Pathway)]

thePlot = ggplot(data_plot,aes(x = abs(NES), y = Pathway, size = -log10(FDR))) + 
  geom_point(shape = 25, fill = "#619CFF") + 
  xlab("Normalized Enrichment") + 
  ylab(NULL) + 
  scale_size_continuous(name = expression(-log[10](FDR)),breaks=1:4) + 
  theme(axis.text.x = element_text(colour = "black")) + 
  theme(axis.text.y = element_text(color = 'black',face=ifelse(data_plot$Pathway == 'Cholesterol\nhomeostasis','bold','plain'))) + 
  theme(legend.position = c(1,0.01),legend.justification = c('right','bottom')) + 
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(plot.margin = margin(5.5, 0, 5.5, 5.5, "pt")) + 
  theme(axis.text.y = element_text(size=10)) +
  theme(legend.title = element_text(size=11),legend.text = element_text(size=10)) + 
  theme(axis.title = element_text(size=12)) + 
  theme(legend.background = element_rect(color='white',fill='white'))

ggsave(filename = sprintf("%s/GSEA_results.pdf",fig_dir), plot = thePlot, width = 8, height = 6)
saveRDS(thePlot,sprintf('%s/GSEA_results.rds',out_dir))