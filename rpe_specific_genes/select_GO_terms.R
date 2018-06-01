library(data.table)
library(stringr)

GSEA_fn = '../processed_data/rpe_specific_genes/GSEA/RPE_hallmark_biocarta_kegg_reactome_GO/gsea_report_for_na_pos_1527185809784.xls'
out_dir = '../processed_data/rpe_specific_genes/select_GO_terms/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}


GSEA = fread(GSEA_fn)
rerun = GSEA[`FDR q-val`==0&str_detect(NAME,'^GO'),NAME]
out_fn = sprintf('%s/rerun.txt',out_dir)
write.table(rerun,out_fn,col.names=FALSE,row.names=FALSE,quote=FALSE)
