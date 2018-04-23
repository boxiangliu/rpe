library(data.table)
library(cowplot)

#-----------#
# Variables #
#-----------#
fig_dir = '../figures/figure4/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

#-----------#
# Functions #
#-----------#
plot_blank = function(){
	ggplot() + geom_blank()
}

make_manhattan_data = function(){
	return(manhattan_data)
}

plot_manhattan = function(manhattan_data){
	return(p)
}

plot_conservation = function(){
	return(p)
}

make_locuszoom_data = function(){
	return(locuszoom_data)
}

plot_locuszoom = function(locuszoom_data){
	return(p)
}

plot_gene_model = function(){
	return(p)
}

#------#
# Main #
#------#

# Manhattan:
blank_plot = plot_blank()

amd_manhattan_data = make_manhattan_data()
amd_manhattan_plot = plot_manhattan(amd_manhattan_data)

myopia_manhattan_data = make_manhattan_data()
myopia_manhattan_plot = plot_manhattan(myopia_manhattan_data)

# Conservation:
conservation = plot_conservation()

# Locuszoom:
eQTL_locuszoom_data = make_locuszoom_data()
eQTL_locuszoom = plot_locuszoom(eQTL_locuszoom_data)

sQTL_locuszoom_data = make_locuszoom_data()
sQTL_locuszoom = plot_locuszoom(sQTL_locuszoom_data)

amd_locuszoom_data = make_locuszoom_data()
amd_locuszoom = plot_locuszoom(amd_locuszoom_data)

myopia_locuszoom_data = make_locuszoom_data()
myopia_locuszoom = plot_locuszoom(myopia_locuszoom_data)

gene_model = plot_gene_model()

# Make Figure 4:
left = plot_grid(amd_manhattan_plot,myopia_manhattan_plot,conservation,nrow=2,align='v',labels = c('A','B','C'))
right = plot_grid(eQTL_locuszoom,sQTL_locuszoom,amd_locuszoom,myopia_locuszoom,gene_model,nrow=5,align='v',labels=c('D','E','F','G',''))
entire = plot_grid(left,right,align = 'h',nrow=1,labels = '')
fig_fn = sprintf('%s/figure4.pdf',fig_dir)
save_plot(fig_fn,entire,base_width=8,base_height=8)