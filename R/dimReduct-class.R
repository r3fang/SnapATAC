#' @importFrom methods setClass setMethod
NULL
#' An S4 class to represent dimentionality reduction object.
#'
#' Class defines a dim.reduct object.
#'
#' @slot dmat a matrix object that contains reduced dimentions
#' @slot sdev variance for each principal conponents
#' @slot iter iterations used for running dimentionality reduction 
#' @slot character a character object indicates the method used for dimentionality reduction
#' @name dim.reduct-class
#' @rdname dim.reduct-class
dim.reduct <- setClass(
  Class = "dim.reduct",
  slots = list(
    dmat = "matrix",
    sdev = "numeric",
	iter = "numeric",
    method = "character"
  )
)

setMethod(
  f = 'show',
  signature = 'dim.reduct',
  definition = function(object) {
    cat(
      'Dimentionality reduction method:', object@method, '\n',
      'Number of dimensions:', ncol(x = object@dmat), '\n'
    )
  }
)

#' subsetting for dim.reduct objects
#'
#' This function takes a dim.reduct object and returns the subset of the object.
#' @param x A dim.reduct object
#' @param i integer; selected rows
#' @param j integer; selected dimentions
#' @param drop character; 
#' @export
setMethod("[", "dim.reduct",
	function(x,i,j, drop="missing"){
		.dmat = x@dmat;
		.sdev = x@sdev;		
		# a single row or column
       if(!missing(i)){
		   if(max(i) > nrow(.dmat)){
			   stop("idx exceeds number of cells");
		   }
		   if(nrow(.dmat) > 0){.dmat <- .dmat[i,,drop=FALSE]}
	   }
	   if(!missing(j)){
		   if(max(j) > ncol(.dmat)){
			   stop("idy exceeds number of dimentions");
		   }
		   if(length(.sdev) > 0){.sdev <- .sdev[j];}	   
 	 	   if(ncol(.dmat) > 0){.dmat <- .dmat[,j,drop=FALSE]}
	   }
	   x@dmat = .dmat;
	   x@sdev = .sdev;
	   return(x);
})


#' @importFrom methods new
newDimReduct <- function () {
	res = new("dim.reduct", 
			  dmat=matrix(0,0,0), 
			  sdev=numeric(),
			  method=character(),
			  iter=numeric()
			  )	
	return(res)
}


isDimReductComplete <- function (obj) {
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "dim.reduct")){
			stop("obj is not a dim.reduct")
		}
	}
	
	if(!((nrow(obj@dmat) > 0) && (length(obj@sdev) > 0))){
		return(FALSE);
	}	
	
	if(ncol(obj@dmat) != length(obj@sdev)){
		return(FALSE);
	}
	
	return(TRUE);
}

weightDimReduct <- function(obj, pca.dims, weight.by.sd=TRUE){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "dim.reduct")){
			stop("obj is not a dim.reduct")
		}
	}
	
	
	if(weight.by.sd){
		data.use = obj@dmat[,pca.dims] %*% diag(sqrt(obj@sdev[pca.dims])) ;
	}else{
		data.use = obj@dmat[,pca.dims];
	}
	return(data.use);
}

#' @importFrom methods is
dimReductDim <- function(obj){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "dim.reduct")){
			stop("obj is not a dim.reduct")
		}
	}
	return(length(obj@sdev));
}
