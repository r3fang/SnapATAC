#' @include utilities.R
globalVariables(names = 'i', package = 'SnapATAC', add = TRUE)
#' Calcualte Jaccard Index Matrix
#'
#' This function takes a snap object as input and calculates the jaccard index matrix 
#' in which each entry Jij equals the intersection over the union between cell i and cell j. 
#' 
#' Calculating jaccard index becomes exponentially time-consuming and also memory 
#' intense with the increase of cell number. 
#' 
#' The most memory and time consuming step of our procedure is the calculation of jaccard index matrix. 
#' This step increases exponentially with the increase of cell number. Here, we calculate a partial jaccard 
#' index matrix against a random subset of reference cells in lieu of the full spectrum of features. 
#' Here we tested if we can perform the calculation of partial jaccard index matrix using a random 
#' subset of max.var reference cells. 
#' 
#' @param obj A Snap obj
#' @param bin.downsample Percentage of bins to be downsampled to [1].
#' @param tmp.folder A non-empty character vector giving the directory name that saves the temp files
#' @param mat A character class that indicates what matrix slot is used to calculate jaccard index c("bmat", "pmat", "gmat")
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted
#' @param seed.use A numeric variable indicates the random seed to use [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(obj=demo.sp, mat="bmat", bin.downsample=1, tmp.folder=tempdir())
#'
#' @return Returns a Snap obj with the jaccard index matrix stored in obj@jmat
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm approxfun density
#' @importFrom methods slot
#' @export
runJaccard <- function(obj, bin.downsample, tmp.folder, mat, max.var, seed.use) {
  UseMethod("runJaccard", obj);
}

#' @export
runJaccard.default <- function(
	obj, 
	tmp.folder,
	bin.downsample=1,
	mat = c("bmat", "pmat", "gmat"),
	max.var = 1000, 
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
	p1 = Matrix::rowMeans(mat.use);

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
	
	# remove the colums in mat.use that have 0 coverage
	col.covs = log(Matrix::colSums(mat.use)+1, 10);
	mat.use = mat.use[,which(col.covs > 0)];
	col.covs = col.covs[which(col.covs > 0)];

	# randomly select a subset of cells as reference 
	if(bin.downsample > 1 | bin.downsample < 0){
		stop("bin.downsample must be between 0 and 1");
	}else{
		if(bin.downsample < 1){
			col.covs.dens <- density(x = col.covs, bw = 'nrd', adjust = 1)
			sampling_prob <- 1 / (approx(x = col.covs.dens$x, y = col.covs.dens$y, xout = col.covs)$y + .Machine$double.eps)
			set.seed(seed.use);
			idy <- sort(sample(x = seq(col.covs), size = bin.downsample*length(col.covs), prob = sampling_prob));
			mat.use = mat.use[,idy];				
		} 
	}
	
	# down-sample for jaccard index matrix
	if(max.var < nrow(obj)){
		max.var = min(max.var, nrow(mat.use));
		row.covs = log(Matrix::rowSums(mat.use)+1,10);		
		row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1)
		sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps)
		set.seed(seed.use);
		idx <- sort(sample(x = seq(row.covs), size = max.var, prob = sampling_prob));
		mat.ref = mat.use[idx,];				
		p2 = p1[idx];		
	}else{
		mat.ref = mat.use;				
		p2 = p1;				
	}
	jmat = calJaccard(mat.use, mat.ref);

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


# Calculate jaccard index between two matrix
# @param X_i first matrix
# @param X_j second matrix
calJaccard <- function(X_i, X_j){
	A = Matrix::tcrossprod(X_i, X_j);
	bi = Matrix::rowSums(X_i);
	bj = Matrix::rowSums(X_j);
	jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
	return(jmat);				
}
