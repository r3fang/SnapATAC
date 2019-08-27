#' @include utilities.R
globalVariables(names = 'i', package = 'SnapATAC', add = TRUE)
#' Normalize Jaccard Index Matrix
#'
#' This function takes a Snap obj as input with jmat slot and normalize 
#' for read depth effect.
#' 
#' In theory, the entries in the jaccard index calculated by calJaccard() should 
#' reflects the true similarity between two cells, however, that is not the case. 
#' We observed that a cell of higher coverage tends to have a higher similarity 
#' with another cell regardless whether these two cells are similar or not.
#' These biases, we termed as “coverage bias” also observed in other studies, 
#' can later result in misleading cell grouping. Therefore, it is cruicial to 
#' normalize the bias. 
#' 
#' @param obj A snap obj
#' @param tmp.folder A non-empty character vector giving the directory name to save temp files.
#' @param method A character class that indicates the normalization method to be used. This must be one of c("residual", "zscore").
#' @param row.center A logical value indicating whether rows of the normalized jaccard inex matrix should be centered by subtracting the layer means (omitting 'NA's).
#' @param row.scale A logical value indicating whether rows of the normalized jaccard index matrix should be scaled by dividing the (centered) layers of 'x' by their standard deviations if 'center' is 'TRUE'.
#' @param high.threshold A numeric class that indicates the max value for normalized jaccard index [5].
#' @param low.threshold A numeric class that indicates the min value for normalized jaccard index [-5].
#' @param do.par A logical variable indicates if to run this in parallel using multiple processors [TRUE].
#' @param ncell.chunk A numeric class that indicates number of cells to process per CPU node
#' @param num.cores A numeric class that indicates the number of cores to use for calculation [1].
#' @param seed.use A numeric class that indicates random seeding number [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runJaccard(obj=demo.sp, mat="bmat", do.par=FALSE, tmp.folder=tempdir());
#' demo.sp = runNormJaccard(obj=demo.sp, do.par=FALSE, tmp.folder=tempdir());
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm
#' @importFrom bigmemory as.big.matrix attach.big.matrix
#' @export
runNormJaccard <- function(obj, tmp.folder, method, row.center, row.scale, low.threshold, high.threshold, do.par, ncell.chunk, num.cores, seed.use){
  UseMethod("runNormJaccard");
}

#' @export
runNormJaccard.default <- function(
	obj, 
	tmp.folder,
	method=c("residual", "zscore"), 
	row.center=TRUE,
	row.scale=TRUE, 
	low.threshold=-5, 
	high.threshold=5, 
	do.par=FALSE,
	ncell.chunk=1000, 
	num.cores=1,
	seed.use=10
){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("'obj' is not a snap obj")
		}
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}

	if(!isJaccardComplete(obj@jmat)){
		stop("jaccard object is not complete, run 'runJaccard' first")
	}else{
		if(isJaccardNorm(obj@jmat)){
			stop("jaccard index matrix has been normalized")
		}		
	}
	
	if(!is.logical(row.center)){
		stop("row.center is not a logical")
	}

	if(!is.logical(row.scale)){
		stop("row.scale is not a logical")
	}
	
	if(low.threshold > high.threshold){
		stop("low.threshold must be smaller than high.threshold");
	}

	if(low.threshold > 0 || high.threshold < 0){
		stop("low.threshold must be smaller than 0 and high.threshold must be greater than 0");
	}
	
	method = match.arg(method);
	jmat = obj@jmat@jmat;
	b1 = obj@jmat@p1;
	b2 = obj@jmat@p2;

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
	
		# step 2) slice the orginal obj into list
		id = seq(nrow(obj));
		id.ls = split(id, ceiling(seq(id)/ncell.chunk));
		
		if(length(id.ls) > 1){
			id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
			# remove the last item of the list
			id.ls = id.ls[-length(id.ls)];
		}	
		
		prefix_tmp = tempfile(pattern = "file", tmpdir = tmp.folder);
		backingfile_tmp <- paste(prefix_tmp, ".bin", sep="");
		descriptorfile_tmp <- paste(prefix_tmp, ".desc", sep="");
	
		x <- as.big.matrix(x = obj@jmat@jmat, 
						   type = "double", 
		                   separated = FALSE, 
						   backingpath=tmp.folder,
		                   backingfile = basename(backingfile_tmp), 
		                   descriptorfile = basename(descriptorfile_tmp)
						   );
	
		cl <- makeCluster(num.cores);
		registerDoParallel(cl);	
		
		nmat <- foreach(i=1:length(id.ls), .verbose=FALSE, .packages="bigmemory", .combine = "rbind") %dopar% {
		    t_mat <- attach.big.matrix(descriptorfile_tmp);
			return(normObservedJmat2(jmat=t_mat[id.ls[[i]],], b1=b1[id.ls[[i]]], b2=b2, method=method));
		}
		
		stopCluster(cl);
		closeAllConnections();
		rm(x);
		file.remove(backingfile_tmp);
		file.remove(descriptorfile_tmp);
		gc();
	}else{
		model.init = trainRegressModel(jmat, b1, b2);
		nmat = normObservedJmat(obj@jmat@jmat, model.init, obj@jmat@p1, obj@jmat@p2, method=method);
	}
	
	if(row.center || row.scale){
		nmat = t(scale(t(nmat), center=row.center, scale=row.scale));
	}
	
	nmat[nmat >= high.threshold] = high.threshold;
	nmat[nmat <= low.threshold]  = low.threshold;

	obj@jmat@jmat = nmat;
	obj@jmat@method = method;
	obj@jmat@norm = TRUE;	
	return(obj);
}

.normOVE <- function(p1, p2){
    pp = tcrossprod(p1, p2);
	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
	ee = pp/(ss - pp)
	return(ee)	
}

trainRegressModel <- function(jmat, b1, b2){
	# remove the diag elements in the jmat
	idx.pairwise = which(jmat == 1, arr.ind=TRUE);

	# calculate the expected jaccard index matrix given the read depth
	emat = .normOVE(b1, b2);

	# estimate the global scaling factor
	scale.factor = mean(jmat / emat);

	# fill the missing value for the diagnoal elements
	jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
	data = data.frame(x=c(emat), y=c(jmat));	
	# 2. polynomial regression
	model <- lm(y ~ x + I(x^2), data);
	return(model);	
}

normObservedJmat <- function(jmat, model, b1, b2, method){
	.normOVE2 <- function(p1, p2){
	    pp = tcrossprod(p1, p2);
		ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
		ee = pp/(ss - pp)
		return(ee)	
	}
	# 1. remove the "1" elements in the jaccard matrix
	idx.pairwise = which(jmat == 1, arr.ind=TRUE);
	emat = .normOVE2(b1, b2);
	scale.factor = mean(jmat / emat);
	jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
	
	# 2. Expansion parameters from subset of cells to all cells
	preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE)
	
	# 3. calculate residuals or zscore
	if(method == "zscore"){
		norm = (c(jmat) - preds$fit) / (preds$se.fit);		
	}else if(method == "residual"){
		norm = c(jmat) -  preds$fit;
	}
	nmat = matrix(norm, nrow(emat), ncol(emat));	
	return(nmat);
}


normObservedJmat2 <- function(jmat, b1, b2, method){
	.normOVE2 <- function(p1, p2){
	    pp = tcrossprod(p1, p2);
		ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
		ee = pp/(ss - pp)
		return(ee)	
	}
	
	# remove the diag elements in the jmat
	idx.pairwise = which(jmat == 1, arr.ind=TRUE);

	# calculate the expected jaccard index matrix given the read depth
	emat = .normOVE(b1, b2);

	# estimate the global scaling factor
	scale.factor = mean(jmat / emat);

	# fill the missing value for the diagnoal elements
	jmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
	data = data.frame(x=c(emat), y=c(jmat));	
	
	# 2. polynomial regression
	model <- lm(y ~ x + I(x^2), data);
	
	# 2. Expansion parameters from subset of cells to all cells
	preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE);
	
	# 3. calculate residuals or zscore
	if(method == "zscore"){
		norm = (c(jmat) - preds$fit) / (preds$se.fit);		
	}else if(method == "residual"){
		norm = c(jmat) -  preds$fit;
	}
	nmat = matrix(norm, nrow(emat), ncol(emat));	
	return(nmat);
}

