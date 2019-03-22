#' @include utilities.R
globalVariables(names = 'i', package = 'SnapATAC', add = TRUE)
#' Normalize Jaccard Index Matrix
#'
#' This function takes a Snap obj as input with jmat slot and normalize 
#' for read depth effect.
#' 
#' In theory, the entry in the jaccard index calculated by calJaccard() should 
#' reflects the true similarity between two cells, however, that is not the case. 
#' We observed that a cell of higher coverage tends to have a higher similarity 
#' with another cell regardless whether these two cells are similar or not.
#' These biases, we termed as “coverage bias” also observed in other studies, 
#' can later result in misleading cell grouping. Therefore, it is cruicial to 
#' normalize the bias. 
#' 
#' Here we propose three different normalization methods.
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
#' @param tmp.folder A non-empty character vector giving the directory name to save temp files
#' @param ncell.chunk A numeric class that indicates number of cells to process per CPU node
#' @param method A character class that indicates the normalization method to be used. This must be one of c("normOVN", "normOVE")
#' @param k A numeric class that indicate number of neibouring cells to use for OVN (used only if method="OVN") (default k = 15) 
#' @param row.center A logical value indicating whether rows of the normalized jaccard inex matrix should be centered by subtracting the layer means (omitting 'NA's)
#' @param row.scale A logical value indicating whether rows of the normalized jaccard index matrix should be scaled by dividing the (centered) layers of 'x' by their standard deviations if 'center' is 'TRUE'.
#' @param high.threshold A numeric class that indicates the max value for normalized jaccard index [5].
#' @param low.threshold A numeric class that indicates the min value for normalized jaccard index [-5].
#' @param num.cores A numeric class that indicates the number of cores to use for calculation [1].
#' @param seed.use A numeric class that indicates random seeding number [10].
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm
#' @importFrom bigmemory as.big.matrix attach.big.matrix
#' @export
#'
runNormJaccard <- function(obj, tmp.folder, ncell.chunk, method, k, row.center, row.scale, low.threshold, high.threshold, num.cores, seed.use){
  UseMethod("runNormJaccard");
}

#' @export
runNormJaccard.default <- function(
	obj, 
	tmp.folder,
	ncell.chunk=1000, 
	method=c("normOVE", "normOVN"), 
	k=15,
	row.center=TRUE,
	row.scale=TRUE, 
	low.threshold=-5, 
	high.threshold=5, 
	num.cores = 1,
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
	
	if(ncell.chunk < 1000){
		stop("ncell.chunk must be larger than 1000")
	}

	if(k < 10 || k > 50){
		stop("k must be in the range between 5 and 50")
	}

	if(low.threshold > high.threshold){
		stop("low.threshold must be smaller than high.threshold");
	}

	if(low.threshold > 0 || high.threshold < 0){
		stop("low.threshold must be smaller than 0 and high.threshold must be greater than 0");
	}
		
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
	
	method = match.arg(method);
	
	message("Epoch: splitting obj into chunks ...");
	# step 2) slice the orginal obj into list
	id = seq(nrow(obj));
	id.ls = split(id, ceiling(seq(id)/ncell.chunk ));
	
	if(length(id.ls) > 1){
		id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
		# remove the last item of the list
		id.ls = id.ls[-length(id.ls)];
	}	
	
	prefix_tmp = tempfile(pattern = "file", tmpdir = tmp.folder);
	backingfile_tmp <- paste(prefix_tmp, ".bin", sep="");
	descriptorfile_tmp <- paste(prefix_tmp, ".desc", sep="");
	
	x <- as.big.matrix(x = obj@jmat@jmat, type = "double", 
	                   separated = FALSE, 
					   backingpath=tmp.folder,
	                   backingfile = basename(backingfile_tmp), 
	                   descriptorfile = basename(descriptorfile_tmp)
					  );

	b1 <- obj@jmat@p1;
	b2 <- obj@jmat@p2;
	
	message("Epoch: scheduling CPUs ...");
	cl <- makeCluster(num.cores);
	registerDoParallel(cl);	

	message("Epoch: normalizing jaccard index for each chunk ...");
	nmat <- foreach(i=1:length(id.ls), .verbose=FALSE, .export="normJaccard", .packages="bigmemory", .combine = "rbind") %dopar% {
	    t_mat <- attach.big.matrix(descriptorfile_tmp);
		return(normJaccard(jmat=t_mat[id.ls[[i]],], b1=b1[id.ls[[i]]], b2=b2, method, k));
	}
	stopCluster(cl);
	closeAllConnections();
	
	message("Epoch: scaling normalized jaccard index ...");
	if(row.center || row.scale){
		nmat = t(scale(t(nmat), center=row.center, scale=row.scale));
	}
	
	message("Epoch: removing values beyond low.threshold and high.threshold ...");
	nmat[nmat >= high.threshold] = high.threshold;
	nmat[nmat <= low.threshold]  = low.threshold;

	message("Epoch: updating jaccard object ...");
	obj@jmat@jmat = nmat;
	obj@jmat@method = method;
	obj@jmat@norm = TRUE;	

	message("Epoch: cleaning up ...");
	rm(x);
	file.remove(backingfile_tmp);
	file.remove(descriptorfile_tmp);
	gc();
	return(obj);
}


