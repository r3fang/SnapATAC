#' Affinity Graph-based Smoothing
#'
#' This function takes a snap object as input with bmat/pmat/gmat/mmat slot 
#' and run magic - an affinity graph-based method - to smooth the signal. 
#'
#' @param obj A snap obj
#' @param input.mat Input matrix c("bmat", "pmat", "gmat", "mmat").
#' @param step.size Number of steps for diffusion [3].
#'
#' @import Matrix
#' @export
runMagic <- function(
	obj,
	input.mat=c("gmat", "pmat", "bmat", "mmat"),
	step.size=3
){
	message("Epoch: checking the inputs ...");
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
	}
	ncell = nrow(obj);
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use = obj@bmat;
		peak.use = obj@feature;
	}else if(input.mat == "pmat"){
		data.use = obj@pmat;
		peak.use = obj@peak;
	}else if(input.mat == "gmat"){
		data.use = obj@gmat;		
		peak.use = GenomicRanges::GRanges();
	}else{
		data.use = obj@mmat;	
		peak.use = GenomicRanges::GRanges();	
	}
	
	if((x=nrow(data.use)) == 0L){
		stop("input matrix is empty")
	}else{
		if((x=nrow(data.use)) != ncell){
			stop("input matrix has wrong number of rows with number of barcodes")
		}
	}
		
	
	# 3. check if the KNN graph exists
	A = obj@graph@mat;
	if((x=nrow(A)) != ncell){
		stop("affnity graph is empty, runKNN first!")
	}
	
	# 4. smooth
	A = A + t(A);
	A = A / Matrix::rowSums(A);
	step.size = 3;
	if(step.size > 1){
	    for(i in 1:step.size){
	        A = A %*% A;
	    }
	}
	data.use.smooth = A %*% data.use;
	slot(obj, input.mat) = data.use.smooth;
	return(obj)
}
