library(data.table)
library(cowplot)
source('utils/gtex_tissue_info.R')

fig_dir = '../figures/figure3/figure3/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

p1 = readRDS('../processed_data/disease_enrichment/mendelian/GEDi_enrichment/GEDi_enrichment.rds')
p1 = p1 + ggtitle('IRD')

p2 = readRDS('../processed_data/disease_enrichment/ld_score_regression/plot_heritability_enrichment/amd_coefficient_pval.rds')
p2 = p2 + ylab(expression(-log[10]*(p-value))) + ggtitle('AMD')

p3 = readRDS('../processed_data/disease_enrichment/ld_score_regression/plot_heritability_enrichment/myopia_coefficient_pval.rds')
p3 = p3 + ylab(expression(-log[10]*(p-value))) + ggtitle('Myopia')

p = plot_grid(p1,p2,p3,labels=c('A','B','C'),nrow=1)
fig_fn = sprintf('%s/figure3.pdf',fig_dir)
save_plot(fig_fn,p,base_height=6,base_width=11)