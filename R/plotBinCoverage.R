#' Plot bin coverage
#'
#' Plot bin coverage distribution
#'
#' @export
plotBinCoverage <- function(obj,...) {
  UseMethod("plotBinCoverage", obj);
}

#' @export
plotBinCoverage.default <- function(obj, ...){	
	# check the input
	if(!(class(obj)=="snap")){stop(paste("Error @plotBarcode: obj is not a snap object!", sep=""))};
	cov = colSums(obj@bmat);
	idy = seq(ncol(colSums(obj@bmat)));	
	cov = cov[which(cov > 0)];	
	cov = log10(cov + 1);
	cov = (cov - mean(cov)) / sd(cov);
	hist(cov, ...);
}

