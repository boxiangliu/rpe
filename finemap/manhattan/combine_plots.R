library(cowplot)

amd_eqtl_fn = '../processed_data/finemap/manhattan/manhattan_amd/manhattan.rds'
amd_sqtl_fn = '../processed_data/finemap/manhattan/manhattan_amd_sqtl/manhattan.rds'
myopia_eqtl_fn = '../processed_data/finemap/manhattan/manhattan_myopia/manhattan.rds'
myopia_sqtl_fn = '../processed_data/finemap/manhattan/manhattan_myopia_sqtl/manhattan.rds'
fig_dir = '../figures/finemap/manhattan/combine_plots/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

amd_eqtl = readRDS(amd_eqtl_fn)
amd_eqtl_p = amd_eqtl[[2]]

amd_sqtl = readRDS(amd_sqtl_fn)
amd_sqtl_p = amd_sqtl[[2]]

p = plot_grid(amd_eqtl_p+ggtitle('eQTL'),amd_sqtl_p+ggtitle('sQTL'),nrow=2,labels = c('A','B'))
fig_fn = sprintf('%s/AMD_manhattan_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)

p = plot_grid(amd_eqtl_p+ggtitle('eQTL'),amd_sqtl_p+ggtitle('sQTL'),nrow=2,labels = c('a','b'))
fig_fn = sprintf('%s/AMD_manhattan_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)


myopia_eqtl = readRDS(myopia_eqtl_fn)
myopia_eqtl_p = myopia_eqtl[[2]]

myopia_sqtl = readRDS(myopia_sqtl_fn)
myopia_sqtl_p = myopia_sqtl[[2]]

p = plot_grid(myopia_eqtl_p+ggtitle('eQTL'),myopia_sqtl_p+ggtitle('sQTL'),nrow=2,labels = c('A','B'))
fig_fn = sprintf('%s/myopia_manhattan_uc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)

p = plot_grid(myopia_eqtl_p+ggtitle('eQTL'),myopia_sqtl_p+ggtitle('sQTL'),nrow=2,labels = c('a','b'))
fig_fn = sprintf('%s/myopia_manhattan_lc.pdf',fig_dir)
save_plot(fig_fn,p,base_height=8,base_width=8)