#' Create Cell-by-Peak Matrix
#'
#' This function takes a snap object and peak list as input and creates cell-by-peak matrix
#' 
#' @param obj A snap object.
#' @param peaks A GRanges object contains the peaks 
#' @param ncell.chunk A numeric class that indicates the number of cells to calculate per processing core [20]. 
#' @param do.par A logical variable indicates if run this using multiple processors [TRUE].
#' @param num.cores Number of processers to use [1].
#' 
#' @return Returns a snap objects contains cell-by-peak matrix
#' 
#' @importFrom plyr count
#' @importFrom GenomicRanges findOverlaps
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @export
createPmat <- function(
	obj,
	peaks,
	ncell.chunk=20,
	do.par=TRUE,
	num.cores=1
){
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap file");			
		}
	}
	
	if(missing(peaks)){
		stop("peaks is missing");
	}else{
		if(!is(peaks, "GRanges")){
			stop("peaks is not a GRanges object");			
		}
	}
	
	if(!is(ncell.chunk, "numeric")){
		stop("ncell.chunk is not numeric object");
	}
	
	if(!is(do.par, "logical")){
		stop("do.par is not logical object");
	}

	if(!is(num.cores, "numeric")){
		stop("num.cores is not logical object");
	}
	
	if(do.par){
	    # input checking for parallel options
		if(num.cores > 1){
	        if (num.cores == 1) {
	          num.cores = 1
	        } else if (num.cores > detectCores()) {
	          num.cores <- detectCores() - 1
	          warning(paste0("num.cores set greater than number of available cores(", parallel::detectCores(), "). Setting num.cores to ", num.cores, "."))
	        }
	      } else if (num.cores != 1) {
	        num.cores <- 1
		}
	}
	
	id = seq(nrow(obj));
	id.ls = split(id, ceiling(seq(id)/ncell.chunk));
	
	barcodes = obj@barcode;
	files = obj@file;
	
	if((x=length(barcodes)) != (y=length(files))){
		stop("obj has different length of barcode and file, recreate the snap object again!")
	}
	
	ind.ls = mclapply(id.ls, function(x){
		reads.gr = extractReads(barcodes[x], files[x], do.par=FALSE, num.cores=1);		
		ov = findOverlaps(reads.gr, peaks);
		idx = match(
			paste(reads.gr$file[queryHits(ov)], reads.gr$barcode[queryHits(ov)]),
			paste(obj@file, obj@barcode)
			);
		idy = subjectHits(ov);		
		ind = data.frame(idx,idy);
		plyr::count(ind, vars = c("idx", "idy"));
	}, mc.cores=num.cores);
	
	ind = Reduce(rbind, ind.ls);
	rm(ind.ls);
	gc();
	
	obj@pmat = 	sparseMatrix(i=ind[,1], 
							 j =ind[,2], 
							 x=ind[,3], 
							 dims=c(length(barcodes), length(peaks))
							 );
	obj@peak = peaks;
	rm(ind);
	gc();
	return(obj);		
}

