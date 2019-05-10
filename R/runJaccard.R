#' @include utilities.R
globalVariables(names = 'i', package = 'SnapATAC', add = TRUE)
#' Calcualte Jaccard Index Matrix
#'
#' This function takes a snap object as input and calculates the jaccard index matrix 
#' in which each entry Jij equals the intersection over the union between cell i and cell j. 
#' 
#' Calculating jaccard index becomes exponentially time-consuming and also memory 
#' intense with the increase of cell number. To solve this problem, we divide
#' the cells into groups and calculated a sub jaccard index matrix, which are later combined
#' to create the full jaccard index matrix.
#' 
#' @param obj A Snap obj
#' @param tmp.folder A non-empty character vector giving the directory name that saves the temp files
#' @param mat A character class that indicates what matrix slot is used to calculate jaccard index c("bmat", "pmat", "gmat")
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted
#' @param ncell.chunk A numeric class that indicates the number of cells to calculate per processing core 
#' @param do.par A logic variable indicates weather to run this in parallel with multiple processors.
#' @param num.cores Number of processors to use.
#' @param seed.use A numeric variable indicates the random seed to use [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(obj=demo.sp, mat="bmat", tmp.folder=tempdir(), do.par=FALSE);
#' 
#' @return Returns a Snap obj with the jaccard index matrix stored in obj@jmat
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm
#' @importFrom methods slot
#' @export

runJaccard <- function(obj, tmp.folder, mat, max.var, ncell.chunk, do.par, num.cores, seed.use) {
  UseMethod("runJaccard", obj);
}

#' @export
runJaccard.default <- function(
	obj, 
	tmp.folder,
	mat = c("bmat", "pmat", "gmat"),
	max.var = 1000, 
	ncell.chunk = 1000, 
	do.par=FALSE,
	num.cores=1,
	seed.use=10
	){
	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		# check if obj is a snap
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
	
	# check if mat is binary;
	mat = match.arg(mat);
	mat.use = methods::slot(obj, mat);
	if((x=nrow(mat.use)) == 0L){
		stop("input matrix is empty");
	}
	
	if((x=max(mat.use)) > 1L){
		stop("input matrix is not a binary matrix, run 'makeBinary' first")	
	}
	
	# check if empty rows exist in bmat;
	if(any(Matrix::rowSums(mat.use) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
    
	max.var = min(max.var, nrow(mat.use));
	# randomly select a subset of cells as reference 
	set.seed(seed.use);
	mat.ref = mat.use[sort(sample(seq(nrow(mat.use)), max.var)),];
	
	if(do.par==TRUE){
	    # input checking for parallel options
		if(num.cores > 1){
	        if (num.cores == 1) {
	          num.cores = 1
	        } else if (num.cores > detectCores()) {
	          num.cores <- detectCores() - 1
	          warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
	        }
	      } else if (num.cores != 1) {
	        num.cores <- 1
		}
		
		# step 2) slice the orginal obj into list
		id = seq(nrow(mat.use));
		id.ls = split(id, ceiling(seq(id)/ncell.chunk));
	
		if(length(id.ls) > 1){
			id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
			# remove the last item of the list
			id.ls = id.ls[-length(id.ls)];
		}		
	
		mat.list = splitBmat(mat.use, id.ls, num.cores, tmp.folder);
	
		cl <- makeCluster(num.cores);
		registerDoParallel(cl);	
	
		jmat <- foreach(i=seq(mat.list), .verbose=FALSE,  .export=c("calJaccardSingle", "calJaccard"), .combine = "rbind") %dopar% {
			calJaccardSingle(mat.list[[i]], mat.ref)
		}
		stopCluster(cl);
		closeAllConnections();
		gc();	
		
		rm(id.ls);
		lapply(mat.list, function(x){
			file.remove(x);
		})
		rm(mat.list);
		gc();
		
	}else{
		jmat = calJaccard(mat.use, mat.ref)
	}

	p1   <- Matrix::rowMeans(mat.use);
	p2   <- Matrix::rowMeans(mat.ref);
	rm(mat.ref);
	rm(mat.use);
	# remove large objects
	obj@jmat@jmat = jmat;
	obj@jmat@p1 = p1;
	obj@jmat@p2 = p2;
	obj@jmat@norm = FALSE;
	obj@jmat@method = character();
	return(obj);
}

