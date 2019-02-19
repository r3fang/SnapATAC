#' Principle Component Anaysis (PCA)
#'
#' This function takes a snap object as input and performs computation of a truncated principal components analysis using irlba
#'
#' @param object A snap object.
#' @param pc.num An integer number of principal component vectors to return [50].
#' @param input.mat A character value indicates what input matrix to use could be c("nmat", "jmat", "cmat", "bmat") [nmat];
#' @param pc.features A numeric index indicates what variables to use [NULL];
#' @param center a logical value indicating whether the variables should be centered to be zero [TRUE]
#' @param scale a logical value indicating whether the variables should be scaled to be standarized [FALSE]
#' @param seed.use a numeric integer value interpreted as an integer [10];
#'
#' @export
runSVD <- function(x, ...) {
  UseMethod("runSVD");
}

runSVD.default <- function(object, pc.num=50, input.mat = c("nmat", "jmat", "cmat", "bmat"), pc.features=NULL, center=TRUE, scale=FALSE, weight.by.var=TRUE, seed.use=10){
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("object is not a snap object")
	}
	ncell = nrow(object)
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "nmat"){
		x = object@nmat;
	}else if(input.mat == "jmat"){
		x = object@jmat;
	}else if(input.mat == "cmat"){
		x = object@cmat;
	}else if(input.mat == "bmat"){
		x = object@bmat;
	}else{
		stop("input.mat does not exist in object")
	}

	if(nrow(x) * ncol(x) == 0){
		stop("input.mat is empty");		
	}

	if(nrow(x) != ncell){
		stop("input.mat has wrong number of cells");
	}
	
	# 3. check pc.features	
	nvar = ncol(x); # number of variables
	if (!is.null(pc.features)) {
		if(any(pc.features > nvar)){stop("pc.features exceeds number of variables in input.mat!")};
		x = x[,pc.features];	
	}
	nvar = ncol(x);

	# 4. check if pc.num is an integer and smaller than cell number;
	if(pc.num%%1!=0){
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
	
	
	# scale the input matrix
	if(center || scale){		
		x <- t(scale(x, center = center, scale = scale));
	}else{
		x <- t(x)
	}
	
	# calculate SVD against the covariance matrix
	set.seed(seed.use);
	S <- irlba(A = x, nv = pc.num, nu = pc.num);
	
	if(weight.by.var){
		object@pmat = S$v %*% diag(S$d);
	}else{
		object@pmat = S$v;
	}
	return(object);
}



