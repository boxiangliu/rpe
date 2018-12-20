library(data.table)
library(stringr)

amd_fn = '../processed_data/finemap/manhattan/2018-06-12_09-41-50_rpe_amd/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'
myopia_fn = '../processed_data/finemap/manhattan/2018-06-06_15-06-09_rpe_23andme/23andme_myopia_prepared_txt_gz_finemap_clpp_status.txt'

amd = fread(amd_fn)
amd[str_detect(feature,'WDR5$'), list(eqtl_file, clpp)] 
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.00177567
# 2:   glucose_eqtls_txt_gz 0.03326801
amd[str_detect(feature,'RLBP1$'), list(eqtl_file, clpp)] 
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.01001784
# 2:   glucose_eqtls_txt_gz 0.00889837

myopia = fread(myopia_fn)
myopia[str_detect(feature,'PDE3A'), list(eqtl_file, clpp)]
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.01353495
# 2:   glucose_eqtls_txt_gz 0.00043300

myopia[str_detect(feature,'PPIL3'), list(eqtl_file, clpp)]
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.00839116
# 2:   glucose_eqtls_txt_gz 0.01015084

myopia[str_detect(feature,'APH1B'), list(eqtl_file, clpp)]
#                 eqtl_file      clpp
# 1: galactose_eqtls_txt_gz 0.0084341
# 2:   glucose_eqtls_txt_gz 0.0622037


myopia[str_detect(feature,'NUCB1'), list(eqtl_file, clpp)]
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.01922592
# 2:   glucose_eqtls_txt_gz 0.01625051

myopia[str_detect(feature,'ETS2'), list(eqtl_file, clpp)]
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.00044152
# 2:   glucose_eqtls_txt_gz 0.01042050

myopia[str_detect(feature,'ENTPD5'), list(eqtl_file, clpp)]
#                 eqtl_file       clpp
# 1: galactose_eqtls_txt_gz 0.00032839
# 2:   glucose_eqtls_txt_gz 0.01340767