gtex_tissue_color_fn = '../data/gtex/gtex_tissue_colors.txt'

read_tissue_color = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_color = gtex_tissue_color$tissue_color_hex
	names(tissue_color) = gtex_tissue_color$tissue_site_detail_id
	tissue_color = c(tissue_color,c(`RPE - glucose` = '#FF0000', `RPE - galactose` = '#00FF00',RPE = '#FF0000'))
	return(tissue_color)
}

read_tissue_abbreviation = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	gtex_tissue_color = gtex_tissue_color[,list(tissue_site_detail_id,tissue_color_hex,tissue_abbreviation)]
	gtex_tissue_color[,tissue_color_hex:=paste0('#',tissue_color_hex)]
	tissue_abbreviation = gtex_tissue_color$tissue_abbreviation
	names(tissue_abbreviation) = gtex_tissue_color$tissue_site_detail_id
	tissue_abbreviation = c(tissue_abbreviation,c(`RPE - glucose` = 'RPE - glucose', `RPE - galactose` = 'RPE - galactose',RPE = 'RPE'))
	return(tissue_abbreviation)
}

get_tissue_to_tissue_id = function(gtex_tissue_color_fn){
	gtex_tissue_color = fread(gtex_tissue_color_fn,colClasses=rep('character',7),sep='\t')
	tissue_to_tissue_id = gtex_tissue_color$tissue_site_detail_id
	names(tissue_to_tissue_id) = gtex_tissue_color$tissue_site_detail
	tissue_to_tissue_id = c(tissue_to_tissue_id,`RPE (glu)`="RPE - glucose",`RPE (gal)`="RPE - galactose",RPE = 'RPE')
	return(tissue_to_tissue_id)
}

tissue_color = read_tissue_color(gtex_tissue_color_fn)
tissue_abbreviation = read_tissue_abbreviation(gtex_tissue_color_fn)
tissue_to_tissue_id = get_tissue_to_tissue_id(gtex_tissue_color_fn)
