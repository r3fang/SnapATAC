#' Feature filtration
#'
#' This function takes a snap obj as input and perform feature selection in the following manner:
#' 1) calculate coverage of each genomic bin/feature;
#' 2) log scale the coverage log10(count + 1);
#' 3) the log-scaled coverage obey approximately a gaussian distribution which is then converted into zscore; 
#' 4) bins with zscore beyond [low.threshold, high.threshold] were filtered;
#' 
#' @param obj A snap obj
#' @param low.threshold Low cutoffs for the parameters (default is -2)
#' @param high.threshold High cutoffs for the parameters (default is 2)
#' @param mat Matrix to filter c("bmat", "pmat")
#'
#' @return Returns a snap obj
#' @importFrom stats sd
#' @importFrom methods slot
#' @export

filterBins <- function(obj, low.threshold, high.threshold, mat) {
  UseMethod("filterBins", obj);
}

#' @export
filterBins.default <- function(obj, low.threshold=-2, high.threshold=2, mat=c("bmat", "pmat")){

	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(class(obj) != "snap"){
			stop("'obj' is not a snap obj")
		};		
	}
	
	mat = match.arg(mat);
	data.use = methods::slot(obj, mat);

	if((x=nrow(data.use)) == 0L){		
		stop("count matrix is empty")
	}

	idy = which(Matrix::colSums(data.use) > 0);
	cov = log(Matrix::colSums(data.use)[idy] + 1, 10);
	zcov = (cov - mean(cov)) / stats::sd(cov);	
	idy2 = which(zcov >= low.threshold & zcov <= high.threshold);
	idy = idy[idy2];
	methods::slot(obj, mat) = data.use[,idy,drop=FALSE];
	return(obj)
}

