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
	}else if(input.mat == "gmat"){
		data.use.ref = obj1@gmat;		
		peak.use.ref = GenomicRanges::GRanges();
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
	
	### ref2
	if(input.mat == "bmat"){
		data.use.qry = obj2@bmat;
		peak.use.qry = obj2@feature;
	}else if(input.mat == "pmat"){
		data.use.qry = obj2@pmat;
		peak.use.qry = obj2@peak;
	}else if(input.mat == "gmat"){
		data.use.qry = obj2@gmat;		
		peak.use.qry = GenomicRanges::GRanges();
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
	
	# subset the matrix for obj2
	message("Epoch: reformatting query dataset input matrix ...");
	idy = unique(queryHits(findOverlaps(peak.use.qry, peak.use.ref)))
	obj2 = obj2[,idy,mat=input.mat]
	# check the dimention of obj1 and obj2
	
	if(!identical(dim(obj1@bmat), dim(obj2@bmat))){
		stop("dimentions of input matrix of obj1 and obj2 do not match")
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










runLDM <- function(
	obj,
	norm=TRUE,
	chr.exclude=NULL,
	bin.filter.cutoff=0.95,
	input.mat=c("bmat", "pmat", "gmat"), 
	num.eigs=20,
	num.landmarks=10000,
	ncell.chunk=10000,
	seed.use=10
){
	nmat.outlier = 0.99;
	message("Epoch: checking the inputs ...");
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
	}
		
	if(!is.logical(norm)){
		stop("norm is not a logical")
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
	
	# exclude the chromsome;
	idy.use = list();
	if(!is.null(chr.exclude)){
		message("Epoch: excluding unwanted chromsomes ...");
		peak.use.copy = peak.use;
		seqlevels(peak.use.copy, force=TRUE) <- setdiff(seqlevels(peak.use.copy), chr.exclude);
		idy.use = append(idy.use, list(which(peak.use$name %in% peak.use.copy$name)));
		rm(peak.use.copy);
	}
	
	bin.cov = log10(Matrix::colSums(data.use)+1);
	cutoff = quantile(bin.cov, bin.filter.cutoff);
	idy.use = append(idy.use, list(which(bin.cov <= cutoff)));
	if(length(idy.use) > 1L){
		idy.use = sort(Reduce(intersect, idy.use));		
	}else{
		idy.use = sort(idy.use[[1]]);
	}
	obj = obj[,idy.use, mat=input.mat];
	
	num.total.cells = nrow(obj);
	if(num.landmarks >= num.total.cells){
		message("Epoch: computing jaccard similarity matrix ...");
		obj = runJaccard2(obj, obj, input.mat=input.mat);
		
		if(norm == TRUE){
			message("Epoch: fitting regression model ...");
			obj = trainRegression(obj)
			message("Epoch: performing normalization ...");
			obj = normJaccard(obj, obj@regModel[1], obj@regModel[2], obj@regModel[3])			
		}else{
			obj@jmat@nmat = obj@jmat@jmat;
		}

		# remove the outliers
		nmat.cutoff = quantile(obj@jmat@nmat, nmat.outlier);
		obj@jmat@nmat[obj@jmat@nmat > nmat.cutoff] = nmat.cutoff;
		
		message("Epoch: computing eigen decomposition ...");
		obj = runEigDecomp(obj, num.eigs);

		# add landmark columns to the metaData
		if(nrow(obj@metaData) == 0){
			obj@metaData$landmark = rep(1, num.total.cells);			
		}else{
			obj@metaData = data.frame(landmark=rep(1, num.total.cells));						
		}
	}else{
		message("Epoch: sampling landmarks ...");
		# calcualte the log-scaled cell coverage
		row.covs = log10(Matrix::rowSums(data.use)+1);
		# estimate the density distribution
		row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1);
		sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps);
		set.seed(seed.use);
		# sample landmarks
		idx.landmark <- sort(sample(x = seq(num.total.cells), size = num.landmarks, prob = sampling_prob));
		idx.other <- setdiff(seq(num.total.cells), idx.landmark);
		obj.landmark = obj[idx.landmark,]

		message("Epoch: computing jaccard similarity matrix for landmarks ...");
		obj.landmark = runJaccard2(obj.landmark, obj.landmark, input.mat=input.mat);
		
		if(norm){
			message("Epoch: fitting regression model ...");
			obj.landmark = trainRegression(obj.landmark)
			message("Epoch: performing normalization ...");
			obj.landmark = normJaccard(obj.landmark, obj.landmark@regModel[1], obj.landmark@regModel[2], obj.landmark@regModel[3]);			
		}else{
			obj.landmark@jmat@nmat = obj.landmark@jmat@jmat;
		}

		# clap the max value for normalize similarity 
		nmat.cutoff = quantile(obj.landmark@jmat@nmat, nmat.outlier);
		obj.landmark@jmat@nmat[obj@jmat@nmat > nmat.cutoff] = nmat.cutoff;
		
		message("Epoch: computing eigen decomposition for landmarks ...");
		obj.landmark = runEigDecomp(obj.landmark, num.eigs);

		idx.other.ls = split(idx.other, ceiling(seq_along(seq(idx.other))/ncell.chunk))
		eig.mat =list()

		message("Epoch: projecting the remaining cells ...");
		
		for(i in seq(idx.other.ls)){
			idx = idx.other.ls[[i]];
			print(paste("Epoch: extending for chunk", i, "...", sep=" "))
			obj.other.i = obj[idx,];
			eig.mat[[i]] = runEigDecompExd(
				obj.other.i, 
				obj.landmark, 
				norm, 
				obj.landmark@regModel[1], 
				obj.landmark@regModel[2], 
				obj.landmark@regModel[3]
			);
			rm(obj.other.i); # free memory
			rm(idx); # free memory
		}
		
		# combine embeddings
		obj@smat@dmat = matrix(0, nr=num.total.cells, nc=num.eigs);
		obj@smat@dmat[idx.landmark,] = obj.landmark@smat@dmat;
		obj@smat@dmat[idx.other,] = as.matrix(do.call(rbind, eig.mat))
		obj@smat@sdev = obj.landmark@smat@sdev

		# add landmark column
		if(nrow(obj@metaData) == 0){
			obj@metaData = data.frame(landmark=rep(0, num.total.cells));
			obj@metaData$landmark[idx.landmark] = 1;
		}else{
			obj@metaData$landmark = rep(0, num.total.cells);
			obj@metaData$landmark[idx.landmark] = 1;
		}
		
		rm(eig.mat);
		rm(obj.landmark);
	}
	if(input.mat == "bmat"){
		obj@bmat = data.use;
		obj@feature = peak.use;
	}else if(input.mat == "pmat"){
		obj@pmat = data.use;
		obj@peak = peak.use;
	}else if(input.mat == "gmat"){
		obj@gmat = data.use;		
	}else{
		stop("input.mat does not exist in obj")
	}
	obj@jmat = newJaccard();
	message("Epoch: Done");
	return(obj);
}
