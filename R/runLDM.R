#' Dimentionality Reduction by Diffusion Maps Algorithm
#'
#' This function takes a snap obj as input and runs diffusion maps 
#' for dimentionality reduction. 
#'
#' Diffusion Maps algorithm, a nonlinear dimensionality reduction technique 
#' that discovers low dimensional manifolds within high-dimensional datasets 
#' by performing harmonic analysis of a random walk constructed over the data 
#' to identify nonlinear collective variables containing the predominance of 
#' the variance in the data. We choose diffusion maps because it is highly robust 
#' to noise and perturbation, making it particuarly suited for analyzing 
#' sparse scATAC-seq dataset.
#' 
#' @param obj A snap obj
#' @param input.mat Input matrix c("bmat", "pmat").
#' @param num.eigs Number of eigenvectors to be computed [20].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runDiffusionMaps(
#'	obj=demo.sp, 
#'	input.mat="bmat", 
#'	num.eigs=20
#' );
#' @import Matrix
#' @export
runDiffusionMaps <- function(
	obj,
	input.mat=c("bmat", "pmat"), 
	num.eigs=20
){
	nmat.outlier = 0.999
	message("Epoch: checking the inputs ...");
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
	}
	
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
	
	rm(data.use); # free memory
	rm(peak.use); # free memory
	
	message("Epoch: computing jaccard similarity matrix ...");
	obj = runJaccard2(obj, obj, input.mat=input.mat);

	message("Epoch: fitting regression model ...");
	obj = trainRegression(obj);

	message("Epoch: performing normalization ...");
	obj = normJaccard(obj, obj@regModel[1], obj@regModel[2], obj@regModel[3]);

	# remove the outliers
	nmat.cutoff = quantile(obj@jmat@nmat, nmat.outlier);
	obj@jmat@nmat[obj@jmat@nmat > nmat.cutoff] = nmat.cutoff;
	
	message("Epoch: computing eigen decomposition ...");
	obj = runEigDecomp(obj, num.eigs);

	obj@smat@method = "DiffusionMaps";
	message("Epoch: Done");
	return(obj);
}


#' Diffusion Maps Extension
#'
#' This function takes two snap objects - one for reference dataset and one for query 
#' dataset and computes the diffusion maps embedding for the query dataset by projecting
#' the query cells into the pre-computed diffusion.
#' 
#' The computational complexity of diffusion maps algorithm exhibits quadratic 
#' growth with the increase of cells, making it infeasible for large-scale datasets. 
#' To overcome this limitation, we apply Nystrom landmark diffusion map algorithm 
#' to efficiently generate the low-dimension embedding for large-scale dataset. 
#' A practical Nystrom landmark diffusion map algorithm project the query dataset 
#' onto the low-dimensional embedding space as learned from the refernce dataset 
#' to create a embedding space for query cells. 
#' 
#' @param obj1 A snap obj for reference dataset
#' @param obj2 A snap obj for query dataset
#' @param input.mat Input matrix c("bmat", "pmat").
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' @import Matrix
#' @export
runDiffusionMapsExtension <- function(
	obj1,
	obj2,
	input.mat=c("bmat", "pmat")
){

	message("Epoch: checking the inputs ...");
	# check if obj1 and obj2 are snap objects
	if(missing(obj1) || missing(obj2)){
		stop("obj1 or obj2 is missing");
	}else{
		if(!is(obj1, "snap") || !is(obj2, "snap")){
			stop("obj1 or obj2 is not a snap obj")
		}		
	}
	
	# check if obj1 has run diffusion maps	
	if((x=nrow(obj1@smat@dmat)) == 0L){
		stop("run diffusion maps 'runDiffusionMaps' to obj1 first!")
	}

	if((x=nrow(obj1@jmat@nmat)) == 0L){
		stop("run diffusion maps 'runDiffusionMaps' to obj1 first!")
	}
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use.ref = obj1@bmat;
		peak.use.ref = obj1@feature;
	}else if(input.mat == "pmat"){
		data.use.ref = obj1@pmat;
		peak.use.ref = obj1@peak;
	}else{
		stop("input.mat does not exist in obj1")
	}
	
	if((x = nrow(data.use.ref)) == 0L){
		stop("input matrix is empty");
	}else{
		if((x = max(data.use.ref)) > 1L){
			stop("input.mat is not a binary matrix")
		}
	}
	
	# check if empty rows exist in bmat;
	if(any(Matrix::rowSums(data.use.ref) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
	
	peak.use.ref.df = as.data.frame(peak.use.ref);
	peak.use.ref$name = paste(peak.use.ref.df[,1], paste(peak.use.ref.df[,2], peak.use.ref.df[,3], sep="-"), sep=":")
	rm(peak.use.ref.df); # free memory 

	### ref2
	if(input.mat == "bmat"){
		data.use.qry = obj2@bmat;
		peak.use.qry = obj2@feature;
	}else if(input.mat == "pmat"){
		data.use.qry = obj2@pmat;
		peak.use.qry = obj2@peak;
	}else{
		stop("input.mat does not exist in obj2")
	}
	
	if((x = nrow(data.use.qry)) == 0L){
		stop("input matrix is empty");
	}else{
		if((x = max(data.use.qry)) > 1L){
			stop("input.mat is not a binary matrix")
		}
	}

	if(any(Matrix::rowSums(data.use.qry) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
	
	peak.use.qry = as.data.frame(peak.use.qry);
	peak.use.qry$name = paste(peak.use.qry[,1], paste(peak.use.qry[,2], peak.use.qry[,3], sep="-"), sep=":")
	rm(peak.use.qry); # free memory 
	
	if(any(peak.use.qry$name != peak.use.ref$name)){
		stop("the column of cell matrix for obj1 and obj2 do not match");
	}
	
	message("Epoch: computing jaccard similarity matrix ...");
	obj2 = runJaccard2(obj2, obj1, input.mat=input.mat);

	message("Epoch: performing normalization ...");
	obj2 = normJaccard(obj2, obj1@regModel[1], obj1@regModel[2], obj1@regModel[3])			

	# remove the outliers
	nmat.cutoff = max(obj1@jmat@nmat);
	obj2@jmat@nmat[obj2@jmat@nmat > nmat.cutoff] = nmat.cutoff;
	
	message("Epoch: projecting query cells to the reference ...");
	obj2 = runEigDecompExd(obj1, obj2);
	
	message("Epoch: Done");
	
	if(input.mat == "bmat"){
		obj2@bmat = data.use.qry
		obj2@feature = peak.use.qry
	}else if(input.mat == "pmat"){
		obj2@peat = data.use.qry
		obj2@peak = peak.use.qry
	}
	return(obj2);
}
