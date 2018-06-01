library(data.table)
library(stringr)

GSEA_fn = '../processed_data/rpe_specific_genes/GSEA/RPE_hallmark_biocarta_kegg_reactome_GO/gsea_report_for_na_pos_1527185809784.xls'
GSEA_rerun_fn = GSEA_fn

read_GSEA = function(fn){
	x = fread(fn,select=c(1,4:11))
	return(x)
}

GSEA = read_GSEA(GSEA_fn)
GSEA = GSEA[str_detect(NAME,'^GO')]

GSEA_rerun = read_GSEA(GSEA_rerun_fn)

replace_rows = function(old,replacement){
	old = copy(old)
	for (i in seq(nrow(replacement))){
		old[NAME == replacement$NAME,] = replacement[i,]
	}
	return(old)
}

GSEA_new = replace_rows(GSEA,GSEA_rerun)
