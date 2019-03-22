#' Normlaize Count Matrix
#'
#' @param obj A snap object
#' @param mat Matrix to be normalized c("gmat", "bmat", "pmat")
#' @param cov A numeric array as cell coverage to normalize the count matrix
#' @param method Method to normalize the matrix c("logRPM", "RPM", "log")
#' importFrom methods as slot
#' @export
scaleCountMatrix <- function(obj, mat, cov, method) {
  UseMethod("scaleCountMatrix", obj);
}

#' @export
scaleCountMatrix.default <- function(
	obj, 
	mat=c("gmat", "bmat", "pmat"), 
	cov, 
	method=c("logRPM", "RPM", "log")
){	
	
	if(!is(obj, "snap")){
		stop("object is not a snap object");
	}

	if(missing(cov)){
		stop("'cov' is not a snap object");
	}
	
	ncell = length(obj@barcode);
	
	if((x = length(cov)) != ncell){
		stop("'cov' has different length from cell number");
	}

	mat = match.arg(mat);
	data.use = methods::slot(obj, mat);
	if((x = nrow(data.use)) == 0){
		stop("selected matrix is empty");
	}
	
	method = match.arg(method);
	if(method == "logRPM"){
		data.use = log((data.use / cov) * 1000000 + 1, 10);
	}else if(method == "RPM"){
		data.use = (data.use) / (cov) * 1000000;		
	}else if(method == "log"){
		data.use = log(data.use+1, 10);		
	}
	
	 methods::slot(obj, mat) = data.use;
	return(obj);
}


