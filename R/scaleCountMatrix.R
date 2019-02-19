#' Normlaize Count Matrix
#'
#' @export

scaleCountMatrix <- function(obj, ...) {
  UseMethod("scaleCountMatrix");
}

scaleCountMatrix.default <- function(obj, mat=c("gmat", "bmat", "pmat"), cov, method=c("RPM", "log")){	
	
	if(class(obj) != "snap"){
		stop("object is not a snap object");
	}

	if(missing(cov)){
		stop("'cov' is not a snap object");
	}
	
	if(length(cov) != length(obj@barcode)){
		stop("'cov' has different length from cell number");
	}

	method = match.arg(method);
	mat = match.arg(mat);
	if(mat == "bmat"){
		if(nrow(obj@bmat) == 0){
			stop("cell-by-bin matrix is empty")			
		}
	}

	if(mat == "pmat"){
		if(nrow(obj@pmat) == 0){
			stop("cell-by-peak matrix is empty")			
		}
	}

	if(mat == "gmat"){
		if(nrow(obj@gmat) == 0){
			stop("cell-by-gene matrix is empty")			
		}
	}

	
	if(mat == "gmat"){
		if(method == "RPM"){
			obj@gmat = (obj@gmat / cov) * 1000000;
		}else{
			obj@gmat = log10(obj@gmat + 1);			
		}
	}else if(mat == "bmat"){
		if(method == "RPM"){
			obj@bmat = (obj@bmat / cov) * 1000000;
		}else{
			obj@bmat = log10(obj@bmat + 1);			
		}
	}else{
		if(method == "RPM"){
			obj@pmat = (obj@pmat / cov) * 1000000;			
		}else{
			obj@pmat = log10(obj@pmat + 1);			
		}		
	}
	return(obj)
}

