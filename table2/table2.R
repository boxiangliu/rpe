library(data.table)
library(foreach)
source('utils/summary_table.R')
source('utils/genome_annotation.R')

glucose_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/glucose/top_intron.txt'
galactose_fn = '../processed_data/sqtl/fastQTL/adjust_pvalue/galactose/top_intron.txt'
out_dir = '../processed_data/table2/table2/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

read_sqtl_result = function(fn){
	sqtl = fread(fn,header=TRUE)
	appendix = foreach(i = seq(nrow(sqtl)),.combine='rbind')%dopar%{
		if (str_detect(sqtl$gene_id[i],',')){
			gene_id = str_split(sqtl$gene_id[i],',')[[1]]
			data.table(gene_id = gene_id,sqtl[i,list(intron,pval,bonf,fdr)])
		} else {
			data.table()
		}
	}
	sqtl = sqtl[!str_detect(gene_id,',')]
	sqtl = rbind(sqtl,appendix)
	return(sqtl)
}

read_gene_annotation = function(){
	gene_annotation = read_gencode(gencode_fn)
	gene_annotation[,type:=NULL]
	setnames(gene_annotation,'gene_type','type')
	return(gene_annotation)
}

select_one_annotation_per_gene = function(sqtl){
	sqtl[,type:=factor(type,levels=c('protein_coding','lincRNA','pseudogene','polymorphic_pseudogene','processed_transcript','sense_intronic','sense_overlapping','antisense','3prime_overlapping_ncrna','unannotated'))]
	sqtl[,rank:=rank(type,ties.method='first'),by='intron']
	sqtl = sqtl[rank==1]
	return(sqtl)
}

gene_annotation = read_gene_annotation()
temp = c(glucose=glucose_fn,galactose=galactose_fn)

for (i in seq_along(temp)){
	sqtl_fn = temp[i]
	sqtl = read_sqtl_result(sqtl_fn)
	sqtl = classify_genes(sqtl, gene_annotation)
	sqtl[is.na(type),type:='unannotated']
	sqtl = select_one_annotation_per_gene(sqtl)

	condensed_sqtl_summary = foreach(FDR = c(1,0.05,0.01,0.001), .combine='rbind')%do%{
		significant_sqtl = sqtl[fdr<=FDR]
		sqtl_class_summary = summarize_genes(significant_sqtl)
		sum(sqtl_class_summary$num)
		condensed_sqtl_summary = condense_summary(sqtl_class_summary, show = c('protein_coding','lincRNA','pseudogene'))
		condensed_sqtl_summary$fdr = FDR
		return(condensed_sqtl_summary)
	}

	reshaped_summary = reshape_summary(condensed_sqtl_summary)
	reshaped_summary = reshaped_summary[condensed_type%in%c('lincRNA','protein_coding'),]
	total = data.frame(condensed_type = 'total',t(colSums(reshaped_summary[,2:ncol(reshaped_summary)])),check.names=FALSE)
	reshaped_summary = rbind(reshaped_summary,total)
	setnames(reshaped_summary,'1','tested')
	for(fdr in c('0.05','0.01','0.001')){
		setnames(reshaped_summary,fdr,'fdr')
		reshaped_summary$pct = calculate_percentage(reshaped_summary$fdr,reshaped_summary$tested)
		reshaped_summary$cbn = sprintf('%s (%s)',reshaped_summary$fdr,reshaped_summary$pct)
		setnames(reshaped_summary,'fdr',fdr)
		setnames(reshaped_summary,'pct',sprintf('%s-pct',fdr))
		setnames(reshaped_summary,'cbn',sprintf('%s-cbn',fdr))
	}

	out_path = sprintf('%s/%s_sqtl_discovery_summary.txt',out_dir,names(temp)[i])
	fwrite(reshaped_summary,out_path,sep='\t')
}
