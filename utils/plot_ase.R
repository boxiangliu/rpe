find_rasqual_fn = function(gene_id,gene_name,dir){
	fn = list.files(dir,pattern=paste0(gene_id,'_',gene_name,'.txt'),recursive=TRUE,full.names=TRUE)
	stopifnot(length(fn)==1)
	return(fn)
}

read_rasqual = function(fn){
	rasqual = fread(fn,header=FALSE)[,c(1,2,3,4,11,12)]
	colnames = c('fid','sid','chr','pos','chisq','pi')
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

extract_ase = function(rasqual){
	ase = rasqual$pi
	pval = rasqual$pval
	zscore = abs(qnorm(pval/2))
	se = abs(ase-0.5)/zscore
	ci = se*1.96
	return(list(ase = ase,se = se,ci = ci))
}

read_ase = function(){
	return(ase)
}

make_ase_dataframe = function(glucose_ase,galactose_ase){
	glucose = data.table(treatment = 'Glucose',ase = glucose_ase$ase, ci = glucose_ase$ci)
	galactose = data.table(treatment = 'Galactose', ase = galactose_ase$ase, ci = galactose_ase$ci)
	ase = rbind(glucose,galactose)
	return(ase)
}

get_ase_plot_data = function(top_tissue_specific_eqtl){
	plot_data=foreach(i = 1:nrow(top_tissue_specific_eqtl),.combine = 'rbind') %do%{
		gene_id = top_tissue_specific_eqtl[i,gene_id]
		gene_name = top_tissue_specific_eqtl[i,gene_name]

		message(gene_id)
		message(gene_name)

		fn = find_rasqual_fn(gene_id,gene_name,glucose_dir)
		rasqual = read_rasqual(fn)
		glucose_ase = extract_ase(rasqual)

		fn = find_rasqual_fn(gene_id,gene_name,galactose_dir)
		rasqual = read_rasqual(fn)
		galactose_ase = extract_ase(rasqual)

		ase = make_ase_dataframe(glucose_ase,galactose_ase)
		ase$gene_name = gene_name
		return(ase)
	}
	plot_data[,gene_name := factor(gene_name,unique(gene_name))]
	plot_data[,treatment:=factor(treatment,levels=c('Glucose','Galactose'))]
	return(plot_data)
}

plot_ase = function(ase){
	ggplot(ase,aes(treatment,ase,ymin = ase-ci,ymax = ase+ci, color = treatment)) + 
		geom_pointrange() + 
		geom_hline(yintercept = 0.5, color = 'red', linetype = 'dashed') + 
		xlab('') + 
		ylab('Allelic ratio (95% CI)') + 
		scale_color_discrete(guide = 'none') +
		facet_grid(.~gene_name) + 
		theme_bw() + 
		theme(strip.background=element_rect(color = 'black', fill = 'white', linetype = 'solid', size = 1))
}

