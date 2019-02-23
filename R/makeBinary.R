#' Convert the count matrix to binary matrix
#'
#' This function takes a snap object as input, then convert count matrix to a binary matrix.
#'
#' It has been observed that some entries in the count matrix have exceedingly high value, which we suspect 
#' are due to alingnment or other type of errors. This function first identifies the top outlier.filter (%) 
#' non-zero values in the count matrix as cutoff and then remove any values above this cutoff as outliers. 
#' 
#' @param object A snap object.
#' @param outlier.filter A numeric class within the range between 0 and 1. The non-zero elememnts in the count matrix are identified annd filtered as outliers (default 1e-3)
#'
#' @return Returns a Snap object with the binary matrix stored in object@bmat
#' 
#' @examples
#' schep.sp
#' schep.sp = cmat2bmat(schep.sp, outlier.filter=1e-3);
#' 
#' @export

makeBinary <- function(object, ...) {
  UseMethod("makeBinary");
}

#' @export
makeBinary.default <- function(object, mat=c("bmat", "pmat", "gmat", outlier.filter=1e-3, )){
	if(class(object) != "snap"){
		stop("object is not a snap object")
	}
	
	if(!is.numeric(outlier.filter) || outlier.filter < 0 || outlier.filter > 1){
		stop("incorrect outlier.filter")		
	}
	
	mat = match.arg(mat);
	if(mat == "bmat"){
		if(nrow(object@bmat) == 0){
			stop("@bmat does not exist")
		}
		
		if(max(object@bmat) == 1){
			stop("@bmat is already binarized")
		}

		x = object@bmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap object
		object@bmat = x;
	}else if(mat == "pmat"){
		if(nrow(object@pmat) == 0){
			stop("@pmat does not exist")
		}
		
		if(max(object@pmat) == 1){
			stop("@pmat is already binarized")
		}

		x = object@pmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap object
		object@pmat = x;
	}else if(mat == "gmat"){
		if(nrow(object@gmat) == 0){
			stop("@gmat does not exist")
		}
		
		if(max(object@gmat) == 1){
			stop("@gmat is already binarized")
		}

		x = object@gmat;
		count = x@x;
		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(count, 1 - outlier.filter));
		# create another binary matrix
		x@x[x@x > count_cutoff] = 0
		x@x[x@x > 0] = 1
		# return snap object
		object@gmat = x;
	}
	return(object);
}

