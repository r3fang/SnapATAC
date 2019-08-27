#' Create Cell-by-Gene Matrix
#'
#' This function takes a snap object and gene list as input and creates cell-by-gene matrix
#' 
#' @param obj A snap object.
#' @param input.mat Input matrix for calculating cell-by-gene matrix.
#' @param genes A GRanges object contains the gene annotation with name as gene name.
#' @param do.par A logical variable indicates if run this using multiple processors [TRUE].
#' @param num.cores Number of processers to use [1].
#' 
#' @return Returns a snap objects contains cell-by-gene matrix stored in "@gmat"
#' 
#' @importFrom plyr count
#' @importFrom GenomicRanges findOverlaps
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @export
createGmatFromMat <- function(obj, input.mat, genes, do.par, num.cores) {
  UseMethod("createGmatFromMat", obj);
}

#' @export
createGmatFromMat.default <- function(
	obj,
	input.mat=c("bmat", "pmat"),
	genes,
	do.par=FALSE,
	num.cores=1
){
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap file");			
		}
	}
	
	if(missing(genes)){
		stop("genes is missing");
	}else{
		if(!is(genes, "GRanges")){
			stop("genes is not a GRanges object");			
		}
		if(is.null(genes$name)){
			stop("genes does not contain gene names");
		}
	}
	
	if(!is(do.par, "logical")){
		stop("do.par is not logical object");
	}

	if(!is(num.cores, "numeric")){
		stop("num.cores is not logical object");
	}
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use = obj@bmat;
		peaks.use = obj@feature;
	}else if(input.mat == "pmat"){
		data.use = obj@pmat;
		peaks.use = obj@peak;
	}else{
		stop("input.mat does not exist in obj")
	}
		
	# check if empty rows exist in bmat;
	if(any(Matrix::rowSums(data.use) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
	num.total.cells = nrow(data.use);
	
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
	
	ovs = as.data.frame(GenomicRanges::findOverlaps(peaks.use, genes));
	ovs.ls = split(ovs, ovs$subjectHits);
	if(do.par){
		count.ls <- parallel::mclapply(ovs.ls, function(idx){
			idx.bins.i = idx$queryHits;
			if(length(idx.bins.i) == 1L){
				count.i = data.use[,idx.bins.i,dropping=TRUE];				
				if(any(count.i > 0)){
					data.frame(i=which(count.i > 0), j=idx$subjectHits[1], val=count.i[count.i > 0])
				}else{
					data.frame()
				}
			}else{
				count.i = Matrix::rowSums(data.use[,idx.bins.i,dropping=TRUE]);				
				if(any(count.i > 0)){
					data.frame(i=which(count.i > 0), j=idx$subjectHits[1], val=count.i[count.i > 0])
				}else{
					data.frame()
				}
			}
		}, mc.cores=num.cores);		
	}else{
		count.ls <- lapply(ovs.ls, function(idx){
			idx.bins.i = idx$queryHits;
			if(length(idx.bins.i) == 1L){
				count.i = data.use[,idx.bins.i,dropping=TRUE];				
				if(any(count.i > 0)){
					data.frame(i=which(count.i > 0), j=idx$subjectHits[1], val=count.i[count.i > 0])
				}else{
					data.frame()
				}
			}else{
				count.i = Matrix::rowSums(data.use[,idx.bins.i,dropping=TRUE]);				
				if(any(count.i > 0)){
					data.frame(i=which(count.i > 0), j=idx$subjectHits[1], val=count.i[count.i > 0])
				}else{
					data.frame()
				}
			}
		});		
	}
	
	count.df = do.call(rbind, count.ls);
	dn = list(obj@barcode, as.character(genes$name))
	obj@gmat = Matrix::sparseMatrix(
		i=count.df[,1], 
		j=count.df[,2], 
		x=count.df[,3], 
		dims=c(nrow(obj), length(genes)),
		dimnames = dn
	);
	rm(count.df);
	rm(dn);
	rm(data.use);
	rm(ovs);
	rm(ovs.ls)
	gc();
	return(obj);		
}
