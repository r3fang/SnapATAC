#' Cell Filteration
#'
#' This function takes a snap object as input and filter nuclei based on given cutoff.
#'
#' @param object snap object
#' @param subset.names Parameters to subset on. Eg, the name of a gene, PC1, a
#' column name in object@@meta.data, etc. Any argument that can be retreived
#' using FetchData
#' @param low.thresholds Low cutoffs for the parameters (default is -Inf)
#' @param high.thresholds High cutoffs for the parameters (default is Inf)
#'
#' @return Returns a snap object containing only the relevant subset of cells
#' 
#' @export
#' 
#' @examples
#' head(x = FetchData(object = pbmc_small, vars.all = 'LTB'))
#' mos_filtered <- filterCells(
#'   object = mos,
#'   subset.names = 'fragment.num',
#'   low.thresholds = 1000
#' )
#' head(x = FetchData(object = pbmc_filtered, vars.all = 'LTB'))
#'
#' @export
filterCells <- function(obj, ...) {
  UseMethod("filterCells");
}

#' @export
filterCells.default <- function(obj, subset.names, low.thresholds, high.thresholds){
	if(class(obj) != "snap"){stop("'obj' is not a snap object")};
	metaData = getMetaData(obj);		
	if(nrow(metaData) == 0){stop("slot metaData is empty")};

    if (missing(x = low.thresholds)){
      low.thresholds <- replicate(n = length(x = subset.names), expr = -Inf);
    }

    if (missing(x = high.thresholds)) {
      high.thresholds <- replicate(n = length(x = subset.names), expr = Inf);
    }

	if(!all(subset.names %in% c('UMI', 'fragment.num', 'mito.ratio', 'umap.ratio', 'dup.ratio', 'pair.ratio'))){
		stop("'subset.names' must be subset of c('UMI', 'fragment.num', 'mito.ratio', 'umap.ratio', 'dup.ratio', 'pair.ratio')");		
	}

    length.check <- sapply(
      X = list(subset.names, low.thresholds, high.thresholds),
      FUN = length
    )
	
    if (length(x = unique(x = length.check)) != 1) {
      stop("'subset.names', 'low.thresholds', and 'high.thresholds' must all have the same length");
    }

	data.subsets <- data.frame(name=subset.names, lw=low.thresholds, hg=high.thresholds);
	
	metaData$UMI = metaData$UQ;
	metaData$fragment.num = metaData$TN;
	metaData$mito.ratio = (metaData$CM)/(metaData$UQ);
	metaData$umap.ratio = (metaData$UM)/(metaData$TN);
	metaData$dup.ratio = 1 - (metaData$UQ)/(metaData$PP);
	metaData$pair.ratio = (metaData$PP)/(metaData$UM);
	
	for(i in seq(nrow(data.subsets))){		
		f = as.character(data.subsets$name[i]);
		names.use <- which(colnames(metaData) %in% f);
		idx = which((metaData[,names.use] >= data.subsets$lw[i]) & (metaData[,names.use] <= data.subsets$hg[i]));
		obj = obj[idx,];
		metaData = metaData[idx,]
	}
	return(obj)
}

