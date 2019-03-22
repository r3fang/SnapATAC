#' Convert the count matrix to binary matrix
#'
#' This function takes a snap obj as input, then convert count matrix to a binary matrix.
#'
#' It has been observed that some entries in the count matrix have exceedingly high value, which we suspect 
#' are due to alingnment or other type of errors. This function first identifies the top outlier.filter (%) 
#' non-zero values in the count matrix as cutoff and then remove any values above this cutoff as outliers. 
#' 
#' @param obj A snap obj.
#' @param mat A character class indicates what matrix slot to use c("bmat", "pmat", "gmat").
#' @param outlier.filter A numeric class within the range between 0 and 1. The non-zero elememnts in the count matrix are identified annd filtered as outliers (default 1e-3)
#' @importFrom stats quantile
#' @return Returns a Snap obj with the binary matrix stored in obj@bmat
#' 
#' @export

makeBinary <- function(obj, mat, outlier.filter) {
  UseMethod("makeBinary", obj);
}

#' @export
makeBinary.default <- function(obj, mat=c("bmat", "pmat", "gmat"), outlier.filter=1e-3){
	if(!is(obj, "snap")){
		stop("obj is not a snap obj")
	}
	
	if(!is.numeric(outlier.filter) || outlier.filter < 0 || outlier.filter > 1){
		stop("incorrect outlier.filter")		
	}
	
	mat = match.arg(mat);
	if(mat == "bmat"){
		if(nrow(obj@bmat) == 0){
			stop("@bmat does not exist")
		}
		
		if(max(obj@bmat) == 1){
			stop("@bmat is already binarized")
		}

		x = obj@bmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap obj
		obj@bmat = x;
	}else if(mat == "pmat"){
		if(nrow(obj@pmat) == 0){
			stop("@pmat does not exist")
		}
		
		if(max(obj@pmat) == 1){
			stop("@pmat is already binarized")
		}

		x = obj@pmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap obj
		obj@pmat = x;
	}else if(mat == "gmat"){
		if(nrow(obj@gmat) == 0){
			stop("@gmat does not exist")
		}
		
		if(max(obj@gmat) == 1){
			stop("@gmat is already binarized")
		}

		x = obj@gmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap obj
		obj@gmat = x;
	}
	return(obj);
}

