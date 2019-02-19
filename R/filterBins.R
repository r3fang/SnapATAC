#' Feature Filteration
#'
#' This function takes a snap object as input and filter feature based on given cutoff.
#'
#' @param object snap object
#' @param low.threshold Low cutoffs for the parameters (default is -2)
#' @param high.threshold High cutoffs for the parameters (default is 2)
#' @param mat input matrix 
#'
#' @return Returns a snap object containing only the relevant subset of bins
#' 
#' @export
#' 

filterBins <- function(object, ...) {
  UseMethod("filterBins", object);
}

#' @export
filterBins.default <- function(object, low.threshold=-2, high.threshold=2, mat=c("bmat", "pmat")){
	if(class(object) != "snap"){
		stop("'object' is not a snap object")
	};
	
	if(mat == "bmat"){
		x = object@bmat;
		if(nrow(x) == 0){
			stop("'object' does not contain cell x bin matrix 'bmat'")
		};
		if(max(x) > 1){
			stop("'bmat' is not a binary matrix, convert it to binary using 'makeBinary' first")
		};		
	}else if(mat == "pmat"){
		x = object@pmat;
		if(nrow(x) == 0){
			stop("'object' does not contain cell x peak matrix 'bmat'")
		};
		if(max(x) > 1){
			stop("'pmat' is not a binary matrix, convert it to binary using 'makeBinary' first")
		};		
	}else{
		stop("'mat' does not exist in object")
	}
	idy = which(Matrix::colSums(x) > 0);
	cov = log(Matrix::colSums(x)[idy] + 1, 10);
	zcov = (cov - mean(cov)) / sd(cov);	
	idy2 = which(zcov >= low.threshold & zcov <= high.threshold);
	idy = idy[idy2];
	if(mat == "bmat"){
		object = object[,idy,"bmat"];		
	}else if(mat == "pmat"){
		object = object[,idy,"pmat"];
	}
	return(object)
}

