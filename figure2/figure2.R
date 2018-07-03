library(cowplot)

manhattan_rds = '../processed_data/diff_expression/plot_manhattan/deseq_manhattan.rds'
GSEA_rds = '../processed_data/diff_expression/plot_GSEA/GSEA_results.rds'
fig_dir = '../figures/figure2/figure2/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

manhattan = readRDS(manhattan_rds)
GSEA = readRDS(GSEA_rds)
blank = ggplot() + geom_blank()

bottom = plot_grid(GSEA,blank,ncol=2,rel_widths = c(1,1),labels = c('b','c'))
entire = plot_grid(manhattan,bottom,nrow=2,rel_heights=c(3,4),labels = c('a',''))
fig_fn = sprintf('%s/figure2_missingC.pdf',fig_dir)

save_plot(fig_fn,entire,base_width=8,base_height=7)