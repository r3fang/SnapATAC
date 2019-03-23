#' Linear Dimentionality Reduction
#'
#' This function takes a snap obj as input and performs computation of linear dimentionality reduction method using irlba
#'
#' @param obj A snap obj.
#' @param pc.num An integer number of dimetions to return [50].
#' @param input.mat A character value indicates what input matrix to use - could be c("nmat", "jmat", "bmat", "pmat") [nmat].
#' @param method A character value indicates what linear dimentionality reduction method to use - could be c("svd", "pca.whiten", "pca") [svd].
#' @param center A logical value indicating whether the variables should be centered to be zero [TRUE].
#' @param scale A logical value indicating whether the variables should be scaled to be standarized [FALSE].
#' @param seed.use A numeric integer value interpreted as an integer [10].
#' @param maxit Maximum number of iterations [1000].
#' @param ... Arguments passed to irlba fuction.
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(obj=demo.sp, tmp.folder=tempdir(), mat="bmat");
#' demo.sp = runNormJaccard(obj=demo.sp, tmp.folder=tempdir());
#' demo.sp = runDimReduct(obj=demo.sp, pc.num=10, input.mat="jmat");
#' 
#' @importFrom irlba irlba prcomp_irlba
#' 
#' @export
runDimReduct <- function(obj, pc.num, input.mat, method, center, scale, seed.use, maxit, ...){
  UseMethod("runDimReduct", obj);
}

#' @export
runDimReduct.default <- function(
	obj, 
	pc.num=50,
	input.mat = c("jmat", "gmat", "bmat", "gmat"), 
	method=c("svd", "pca.whiten", "pca"), 
	center=TRUE, 
	scale=FALSE, 
	seed.use=10,
	maxit=1000,
	...
){
	# 1. check if obj is a snap;
	if(!is(obj, "snap")){
		stop("obj is not a snap obj")
	}
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "jmat"){
		x = obj@jmat@jmat;
	}else if(input.mat == "bmat"){
		x = obj@bmat;
	}else if(input.mat == "pmat"){
		x = obj@pmat;
	}else if(input.mat == "gmat"){
		x = obj@gmat;
	}else{
		stop("input.mat does not exist in obj")
	}
	
	if(nrow(x) * ncol(x) == 0){
		stop("input.mat is empty");		
	}
	
	ncell = nrow(obj);
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
		set.seed(seed.use)
		S <- irlba(A = x.norm, nv = pc.num, nu = pc.num, maxit=maxit, scale.=FALSE, center=FALSE, ...);
		obj@smat@dmat = S$v;
		obj@smat@sdev = S$d;
		obj@smat@iter = maxit;
		obj@smat@method = "svd";
	}else if(method == "pca.whiten"){
		# calculate covariance matrix
		V <- x %*% t(x) / ncell;
		# calculate SVD against the covariance matrix
		set.seed(seed.use);
		S <- irlba(A = V, nv = pc.num, nu = pc.num);
		# PCA whitening
		D <- diag(c(1/sqrt(S$d)))
		K <- D %*% t(S$u)
		obj@smat@dmat = t(K %*% x);
		obj@smat@sdev = S$d;
		obj@smat@iter = maxit;
		obj@smat@method = "pca.whiten";		
	}else{
		set.seed(seed.use);
		S <- prcomp_irlba(x.norm, n=pc.num, scale.=FALSE, center=FALSE, maxit=maxit, ...);
		obj@smat@dmat = S$x;
		obj@smat@sdev = S$d;
		obj@smat@iter = maxit;
		obj@smat@method = "pca";
	}
	return(obj);
}



