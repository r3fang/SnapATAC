#' @importFrom methods setClass setMethod
NULL
#' An S4 class jaccard to represent a jaccard object.
#'
#' Class jaccard defines a jaccard object.
#'
#' @slot jmat a matrix object that contains the jaccard index similarity matrix
#' @slot nmat a matrix object that contains the normalized similarity matrix
#' @slot p1 an array of numeric values indicates the coverage for rows
#' @slot p2 an array of numeric values indicates the coverage for columns
#' @slot norm a logical variable indicates whether jaccard index matrix has been normalized
#' @slot input.mat a character variable indicates what matrix is used as input (i.e. "bmat", "pmat", "gmat")
#' @examples
#' data(demo.sp);
#' demo.sub.sp = demo.sp[1:5,]
#' @name jaccard-class
#' @rdname jaccard-class
#' @importFrom methods setClassUnion
methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
jaccard <- setClass(
  Class = "jaccard",
  slots = list(
    jmat = "MatrixOrmatrix",
    nmat = "MatrixOrmatrix",
    p1 = "numeric",
	p2 = "numeric",
	norm = "logical",
	input.mat="character"
  )
)

setMethod(
  f = 'show',
  signature = 'jaccard',
  definition = function(object) {
    cat(
      'Number of cells:', nrow(object@jmat), '\n',
      'Number of dims: ', ncol(object@jmat), '\n',
	  'Input matrix:', object@input.mat, "\n",
      'Normalized: ', object@norm, '\n'
    )
  }
)

#' subsetting for jaccard objects
#'
#' This function takes a jaccard object and returns the subset of jaccard object.
#' @param x A jaccard object
#' @param i selected rows
#' @param j selected columns
#' @param drop drop unused levels
#' @examples
#' data(demo.sp);
#' demo.sub.sp = demo.sp[1:5,]
#' @export
setMethod(
	f="[", 
    signature = 'jaccard',
	function(x,i,j, drop="missing"){
		.jmat = x@jmat;
		.nmat = x@nmat;	
		.p1 = x@p1;		
		.p2 = x@p2;			
		# a single row or column
       if(!missing(i)){
		   if(max(i) > nrow(.jmat)){
			   stop("idx exceeds number of cells");
		   }
		   if(nrow(.jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
		   if(nrow(.nmat) > 0){.nmat <- .nmat[i,,drop=FALSE]}
		   if(length(.p1) > 0){.p1 <- .p1[i,drop=FALSE]}
	   }
	   if(!missing(j)){
		   if(max(j) > ncol(.jmat)){
			   stop("idy exceeds number of dimentions");
		   }
 	 	   if(ncol(.jmat) > 0){.jmat <- .jmat[,j,drop=FALSE]}
		   if(ncol(.nmat) > 0){.nmat <- .nmat[,j,drop=FALSE]}
		   if(length(.p2) > 0){.p2 <- .p2[j,drop=FALSE]}
	   }
	   x@jmat = .jmat;
	   x@p1 = .p1;
	   x@p2 = .p2;
	   x@nmat = .nmat;
	   return(x);
})

#' @importFrom methods new
newJaccard <- function () {
	res = new("jaccard", 
			  jmat=matrix(0,0,0), 
			  p1=numeric(),
			  p2=numeric(),
			  norm=FALSE,
			  nmat=matrix(0,0,0),
			  input.mat="None"
			  )	
}

isJaccardComplete <- function (obj) {
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "jaccard")){
			stop("obj is not a jaccard object")
		}
	}
	if(!((nrow(obj@jmat) > 0) && (length(obj@p1) > 0) && (length(obj@p2) > 0))){
		return(FALSE);
	}	
	return(TRUE);
}


isJaccardNorm <- function (obj) {
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "jaccard")){
			stop("obj is not a jaccard object")
		}
	}
	return(obj@norm)
}

