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
#' @param mat A character class that indicates what matrix slot is used to calculate jaccard index c("bmat", "pmat", "gmat")
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted [2000]
#' @param ncell.chunk A numeric class that indicates the number of cells to calculate per processing core [1000]
#' @param seed.use A numeric variable indicates the random seed to use [10]
#' @return Returns a Snap obj with the jaccard index matrix stored in obj@jmat
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm
#' @importFrom methods slot
#' @export

runJaccard <- function(obj, mat, max.var, ncell.chunk, seed.use) {
  UseMethod("runJaccard", obj);
}

#' @export
runJaccard.default <- function(
	obj, 
	mat = c("bmat", "pmat", "gmat"),
	max.var = 2000, 
	ncell.chunk = 1000, 
	seed.use = 10
	){
	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		# check if obj is a snap
		if(class(obj) != "snap"){
			stop("'obj' is not a snap obj")
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
    
	num.cores <- 1
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
    num.cores <- 1
	# currently this dunction does not support parallele computing
	#num.cores = 1;
	# max.var can not exceed number of cells
	max.var = min(max.var, nrow(mat.use));

	# randomly select a subset of cells as reference 
	set.seed(seed.use);
	mat.ref = mat.use[sort(sample(seq(nrow(mat.use)), max.var)),]
	
	# step 2) slice the orginal obj into list
	id = seq(nrow(mat.use));
	id.ls = split(id, ceiling(seq_along(id)/ncell.chunk));
	
	message("Epoch: splitting obj into chunks ...");
	if(length(id.ls) > 1){
		id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
		# remove the last item of the list
		id.ls = id.ls[-length(id.ls)];
	}	
	
	message("Epoch: calculating jaccard index for each chunk ...");	
	## Save the split data in an environment
	my.env <- new.env();	
	for(i in seq(id.ls)){
		mat_i = mat.use[id.ls[[i]],];		
		assign(paste("mat", i, sep=""), mat_i, envir=my.env)		
	}
	my.env$mat0 = mat.ref;
	
	## Function that takes indexes, then extracts the data from the environment
	applyMyFun <- function(idx, env) {
	    eval(parse(text = paste0("mat_i <- env$mat", idx)))
		calJaccard(mat_i, env$mat0)
	}
	
	message("Epoch: scheduling CPUs ...");
	cl <- makeCluster(num.cores);
	registerDoParallel(cl);	
	
	message("Epoch: calculating jaccard index for each chunk ...");
	jmat <- foreach(i=seq(id.ls), .verbose=FALSE,  .export="calJaccard", .combine = "rbind") %dopar% {
		applyMyFun(i, my.env);
	}
	stopCluster(cl);
	closeAllConnections();
	
	p1   <- Matrix::rowMeans(mat.use);
	p2   <- Matrix::rowMeans(mat.ref);

	# remove large objects
	message("Epoch: cleaning up ...");
	rm(my.env);
	rm(mat.ref);
	rm(id.ls);
	gc();
	obj@jmat@jmat = jmat;
	obj@jmat@p1 = p1;
	obj@jmat@p2 = p2;
	obj@jmat@norm = FALSE;
	obj@jmat@method = character();
	return(obj);
}






