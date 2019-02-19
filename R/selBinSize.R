#' checkBinSize
#'
#' This function takes a snap object as input and check if the current bin size is ok
#'
#' @param obj character. a snap object
#' @export
checkBinSize <- function(obj1, ...) {
  UseMethod("checkBinSize", obj1);
}

#' @export
checkBinSize.default <- function(obj1, obj2=NULL){
	if(!is.snap(obj1)){
		stop("'obj1' is not a snap object")
	}

	if(nrow(obj1@bmat) == 0){
		stop("@bmat cell-by-bin matrix is empty")		
	}
	
	# if obj2 is given
	if(!is.null(obj2)){
		# check if obj2 is a snap object
		if(!is.snap(obj2)){
			stop("'obj2' is not a snap object")
		}
		
		# check if obj2 has the same features
		if(any(obj1@feature$name != obj2@feature$name)){
			stop("'obj1' and 'obj2' have different features")			
		}
		obj1 = makeBinary(obj1);						
		obj2 = makeBinary(obj2);						
		cov1 = log(Matrix::colSums(obj1@bmat) + 1, 10);
		cov2 = log(Matrix::colSums(obj2@bmat) + 1, 10);	
		
	}else{
		obj1 = makeBinary(obj1);						
		ncell = length(obj1@barcode);
		idx1 = sort(sample(seq(ncell), ncell/2));
		idx2 = setdiff(seq(ncell), idx1);
		cov1 = log(Matrix::colSums(obj1@bmat[idx1,]) + 1, 10);
		cov2 = log(Matrix::colSums(obj1@bmat[idx2,]) + 1, 10);	
	}		
	return(cor(cov1, cov2));
}

