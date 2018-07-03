library(manhattan)

plot_manhattan = function(data,label=NULL, cutoff = CUTOFF){
	dummy=data.table(condition=c('Glucose','Galactose'),y=c(cutoff,cutoff))
	p=manhattan(data,build='hg19')+
		facet_grid(condition~.,scale='free_y')+
		scale_y_sqrt(breaks = c(0.01,seq(0.1,0.5,0.1)))+
		geom_text_repel(aes(label=label),color='black')+
		geom_hline(data=dummy,aes(yintercept=y),color='red',linetype=2)+
		ylab(paste('Colocalization probability')) +
		theme_bw() + 
		theme(panel.grid = element_blank(), panel.border = element_rect(size = 1), strip.background = element_rect(color = 'black',fill = 'white', size = 1)) 
	return(p)
}