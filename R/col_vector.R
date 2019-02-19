colVector <- function(n=60, type=c("qual", "div", "seq")){
	library(RColorBrewer);
	type <- match.arg(type);
	set.seed(10)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == type,];
	set.seed(10)
	col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));
	return(col_vector)
} 

