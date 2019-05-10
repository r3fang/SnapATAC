#' Jaccard Distance Analysis
#'
#' This function takes a Snap obj as input with bmat/pmat slot and run 
#' Jaccard distance based analysis (JDA).
#' 
#' runJDA performs the following steps:
#' 1) runJaccard - calculate pair-wise jaccard index
#' 2) normJaccard - adjusts read depth and other artifacts in observed jaccard matrix
#' 3) SVD - run SVD on the normalized jaccard matrix
#' 
#' In theory, the entry in the jaccard index calculated by calJaccard() should 
#' reflects the true similarity between two cells, however, that is not the case. 
#' We observed that a cell of higher coverage tends to have a higher similarity 
#' with another cell regardless whether these two cells are similar or not.
#' These biases, we termed as “coverage bias” also observed in other studies, 
#' can later result in misleading cell grouping. Therefore, it is cruicial to 
#' normalize the bias. 
#' 
#' #' Here we propose three different normalization methods.
#'	
#' 1. Observe Over Neighbours (OVN)
#' 2. Observe Over Expected   (OVE)
#'
#' 1. Observe Over Neighbours (OVN)
#' For each pair of cells i and j, we first identify their neibours Ni and Nj based
#' on the coverage. Second, we calcualted the pair-wise jaccard index between Ni 
#' Nj as Eij. The average mean(Eij) is considered to be the expected jaccard index 
#' between cell i and j. The normalized jaccard index between i and j is Oij - mean(Eij).
#'
#' 2. Observe Over Expected (OVE)
#' Alternatively, one can calcualte the theoritical expected jaccard index between 
#' cell i and j if the "1"s are completely random. However, we found the theoretically
#' estimated expected jaccard index is usually lower than that estimated from the neibours.
#' This is also expected because the accessible sites are not completely random. For instance,
#' we previously found the promoters of housekeeping genes are more often occupied between
#' different cells. To solve this under-estimate problem for OVE, we found that OVE is highly
#' correlated with OVN (r>0.99), therefore, we can estimate a scale factor alpha to scale 
#' OVE to similar level with OVN. The scale factor is estimated from the data.
#' 
#' @param obj A snap obj
#' @param input.mat Input matrix to be used for LSA c("bmat", "pmat", "gmat").
#' @param bin.cov.zscore.lower Bin coverage is coverted to zscore and bins with zscore lower than bin.cov.zscore.lower will be filtered
#' @param bin.cov.zscore.upper Bin coverage is coverted to zscore and bins with zscore higher than bin.cov.zscore.upper will be filtered
#' @param pc.num An integer number of dimetions to return [50].
#' @param norm.method A character class that indicates the normalization method to be used. This must be one of c("normOVE", "normOVN")
#' @param tmp.folder A non-empty character vector giving the directory name that saves the temp files
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted
#' @param row.center A logical value indicating whether rows of the normalized jaccard inex matrix should be centered by subtracting the layer means (omitting 'NA's)
#' @param row.scale A logical value indicating whether rows of the normalized jaccard index matrix should be scaled by dividing the (centered) layers of 'x' by their standard deviations if 'center' is 'TRUE'.
#' @param high.threshold A numeric class that indicates the max value for normalized jaccard index [5].
#' @param low.threshold A numeric class that indicates the min value for normalized jaccard index [-5].
#' @param do.par A logic variable indicates weather to run this in parallel with multiple processors.
#' @param num.cores Number of processors to use.
#' @param ncell.chunk A numeric class that indicates the number of cells to calculate per processing core 
#' @param seed.use A numeric variable indicates the random seed to use [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runJDA(
#'	obj=demo.sp, 
#'	input.mat="bmat", 
#'  bin.cov.zscore.lower=-2,
#'  bin.cov.zscore.upper=2,
#'	pc.num=50,
#'	norm.method="normOVE",
#'	tmp.folder=tempdir(),
#'	max.var=2000,
#'	do.par=TRUE,
#'	ncell.chunk=1000,
#'	num.cores=5,
#'	seed.use=10
#'	);
#' 
#' @import Matrix
#' @export

runJDA <- function(
	obj,
	input.mat = c("bmat", "pmat"),
	bin.cov.zscore.lower=-2,
	bin.cov.zscore.upper=2,
	pc.num=50,
	norm.method=c("normOVE", "normOVN"),
	tmp.folder,
	max.var=2000,
	row.center=TRUE,
	row.scale=TRUE, 
	low.threshold=-5, 
	high.threshold=5,
	do.par=TRUE,
	ncell.chunk=1000,
	num.cores = 1,
	seed.use=10
){
	message("Epoch: checking the inputs ...");
	if(missing(obj)){
		stop("obj is missing");
	}else{
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
	
	
	if(!is.logical(row.center)){
		stop("row.center is not a logical")
	}

	if(!is.logical(row.scale)){
		stop("row.scale is not a logical")
	}
	
	if(ncell.chunk < 1000){
		stop("ncell.chunk must be larger than 1000")
	}
		
	if(low.threshold > high.threshold){
		stop("low.threshold must be smaller than high.threshold");
	}

	if(low.threshold > 0 || high.threshold < 0){
		stop("low.threshold must be smaller than 0 and high.threshold must be greater than 0");
	}

	if(bin.cov.zscore.lower > bin.cov.zscore.upper){
		stop("bin.cov.zscore.lower must be smaller than bin.cov.zscore.upper");
	}
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use = obj@bmat;
	}else if(input.mat == "pmat"){
		data.use = obj@pmat;
	}else{
		stop("input.mat does not exist in obj")
	}
	
	if((x = nrow(data.use)) == 0L){
		stop("input matrix is empty");
	}else{
		if((x = max(data.use)) > 1L){
			stop("input.mat is not a binary matrix")
		}
	}
	# check if empty rows exist in bmat;
	if(any(Matrix::rowSums(data.use) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
	rm(data.use);

	message("Epoch: filtering bins ..");
	obj = filterBins(
		obj=obj, 
		low.threshold=bin.cov.zscore.lower,
		high.threshold=bin.cov.zscore.upper
	);
		
	message("Epoch: running jaccard index matrix ...");
	obj = runJaccard(
		obj=obj, 
		tmp.folder=tmp.folder, 
		mat=input.mat, 
		max.var=max.var, 
		do.par=FALSE,
		seed.use=seed.use
	);	
	
	message("Epoch: normalizing jaccard index matrix ...");
	obj = runNormJaccard(
		obj=obj, 
		tmp.folder=tmp.folder,
		do.par=do.par,
		ncell.chunk=ncell.chunk, 
		method=norm.method,
		k=15,
		row.center=row.center,
		row.scale=row.scale, 
		low.threshold=low.threshold, 
		high.threshold=high.threshold, 
		num.cores =num.cores,
		seed.use=seed.use
	);
	
	message("Epoch: running dimentionality reduction ...");
	obj = runDimReduct(
		obj=obj, 
		pc.num=pc.num,
		input.mat = "jmat",
		method="svd",
		center=TRUE, 
		scale=FALSE, 
		seed.use=seed.use
	);	
	obj@jmat = newJaccard();
	message("Epoch: Done");
	return(obj);
}

