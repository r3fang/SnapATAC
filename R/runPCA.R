#' Linear Dimentionality Reduction
#'
#' This function takes a snap object as input and performs computation of linear dimentionality reduction method using irlba
#'
#' @param object A snap object.
#' @param pc.num An integer number of dimetions to return [50].
#' @param input.mat A character value indicates what input matrix to use - could be c("nmat", "jmat", "bmat", "pmat") [nmat];
#' @param method A character value indicates what linear dimentionality reduction method to use - could be c("svd", "pca.whiten", "pca") [svd]
#' @param center A logical value indicating whether the variables should be centered to be zero [TRUE]
#' @param scale A logical value indicating whether the variables should be scaled to be standarized [FALSE]
#' @param weight.by.sd A logical value indicating whether to scale the vector using standard deviation [FALSE];
#' @param seed.use a numeric integer value interpreted as an integer [10];
#'
#' @export
runPCA <- function(object, ...) {
  UseMethod("runPCA", object);
}

#' @export
runPCA.default <- function(
	object, 
	pc.num=50,
	input.mat = c("nmat", "jmat", "cmat", "bmat"), 
	method=c("svd", "pca.whiten", "pca"), 
	weight.by.sd = TRUE,
	center=TRUE, 
	scale=FALSE, 
	seed.use=10
){
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("object is not a snap object")
	}
	ncell = length(object@barcode)

	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "nmat"){
		x = object@nmat;
	}else if(input.mat == "jmat"){
		x = object@jmat;
	}else if(input.mat == "bmat"){
		x = object@bmat;
	}else if(input.mat == "pmat"){
		x = object@pmat;
	}else{
		stop("input.mat does not exist in object")
	}
	
	if(nrow(x) * ncol(x) == 0){
		stop("input.mat is empty");		
	}

	if(nrow(x) != ncell){
		stop("input.mat has wrong number of cells");
	}
	
	method = match.arg(method);
	nvar = ncol(x); # number of variables

	# 4. check if pc.num is an integer and smaller than cell number;	
	if (is.numeric(pc.num)) {
	    pc.num <- as.integer(pc.num)
	  } else{
  		stop("pc.num must be an integer")	  	
	}

    if (pc.num > min(ncell, nvar)) {
        message("'n.comp' is too large: reset to ", min(ncell, nvar))
        pc.num <- min(nvar, ncell)
    }
	
	# 5. check if scale and center is logical variable
	if(!is.logical(center)){
		stop("center must be logical variable TRUE or FALSEs")
	}
	if(!is.logical(scale)){
		stop("scale must be logical variable TRUE or FALSEs")
	}
	
	x <- t(x)
	
    a <- names(as.list(match.call()))
    ans <- list(scale=scale)
	 
    if (!is.matrix(x)) x <- as.matrix(x)
    args <- list(A=x, nv=pc.num)
    if (is.logical(center))
    {
      if (center) args$center <- colMeans(x)
    } else{
    	stop("center must be logical")
    }
    if (is.logical(scale))
    {
        if (is.numeric(args$center))
        {
          f <- function(i) sqrt(sum((x[, i] - args$center[i]) ^ 2) / (nrow(x) - 1L))
          scale. <- vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE)
          if (ans$scale) ans$totalvar <- ncol(x)
          else ans$totalvar <- sum(scale. ^ 2)
        } else
        {
          if (ans$scale)
          {
            scale. <- apply(x, 2L, function(v) sqrt(sum(v ^ 2) / max(1, length(v) - 1L)))
            f <- function(i) sqrt(sum((x[, i] / scale.[i]) ^ 2) / (nrow(x) - 1L))
            ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE) ^ 2)
          } else
          {
            f <- function(i) sum(x[, i] ^ 2) / (nrow(x) - 1L)
            ans$totalvar <- sum(vapply(seq(ncol(x)), f, pi, USE.NAMES=FALSE))
          }
        }
        if (ans$scale) args$scale <- scale.
    } else{
    	stop("scale must be logical")
    }

	if(center){
		x.norm = sweep(args$A, 2, args$center, FUN=`-`)		
	}else{
		x.norm = args$A;		
	}
	if(scale){
		x.norm = sweep(x.norm, 2, args$scale, FUN=`/`)		
	}else{
		x.norm = x.norm
	}
	
	if(method == "svd"){
		S <- irlba(A = x.norm, nv = pc.num, nu = pc.num);		
		if(weight.by.sd){
			object@smat = S$v %*% diag(sqrt(S$d));
		}else{
			object@smat = S$v;					
		}
		object@d = S$d;
	}else if(method == "pca.whiten"){
		# calculate covariance matrix
		V <- fastMatMult(x, t(x))/ncell;
		# calculate SVD against the covariance matrix
		set.seed(seed.use);
		S <- irlba(A = V, nv = pc.num, nu = pc.num);
		# PCA whitening
		D <- diag(c(1/sqrt(S$d)))
		K <- D %*% t(S$u)
		object@smat = t(K %*% x);		
		object@d = S$d;
	}else{
		S <- prcomp_irlba(x.norm, n=pc.num, scale=FALSE, center=FALSE);
		object@smat = S$x;		
		object@d = S$d;
	}
	return(object);
}



