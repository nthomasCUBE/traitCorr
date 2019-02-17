library(d3heatmap)
library(scales)
library(xlsx)
library(gplots)

options(stringsAsFactors=FALSE)

add_transcriptomics=function(){
	print(c("INFO|add_transcriptomcis"))
}

add_trait_information=function(){
	print(c("INFO|add_trait_information"))
}

calc_cmp_transcriptomics_traits=function(v){
	print(c("INFO|transcriptomics VS traits"))
	print(dim(v$transcriptomics))
	print(dim(v$trait))
}


