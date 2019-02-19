#' Calcualte Jaccard Index Matrix
#'
#' To conquar the challenge of sparsity, we first converted the sparse binary matrix 
#' into a jaccard index matrix in which each entry Jij equals the intersection over the 
#' union between cell i and cell j. The end result is a symmetric matrix with each
#' element representing the similarity between two cells. 
#' 
#' Calculating jaccard index becomes exponentially time-consuming and also memory 
#' intense with the increase of cell number. To solve this problem, we first divided 
#' the cells into groups and calculated a sub jaccard index matrix, which are combined
#' to create the full jaccard index matrix.
#' 
#' For instance, given that there are 50,000 cells in total, we first split the cells into 
#' 10 chunks with each chunk containing 5,000 cells. Then we calculated the pairwise sub 
#' jaccard index matrix for each chunk. Finally, we created the entire jaccard index matrix 
#' by combining all sub jaccard matrices. 
#' 
#' To further speed up the calculation and avoid of storing a huge unsparse matrix, 
#' we allow for calculating a partial jaccard matrix.
#'
#' @param object A Snap object
#' @param cell.pcore A numeric class that indicates the number of cells to calculate per processing core [1000]
#' @param ncore A numeric class that indicates the number of cores to use for calculation [1]
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted [5000]
#' @param norm A logic variable indicates the whether to normalize the jaccard index matrix [TRUE]
#' @param seed.use A numeric variable indicates the random seed to use [10]
#'
#' @return Returns a Snap object with the jaccard index matrix stored in object@jmat
#' @import doSNOW
#' @export

calJaccard <- function(object, ...) {
  UseMethod("calJaccard", object);
}

#' @export
calJaccard <- function(
	object, 
	mat = c("bmat", "pmat", "gmat"),
	ncell.chunk=5000, 
	max.var=5000, 
	seed.use=10, 
	norm.method=c("normOVE", "normOVN"), 
	k=15,
	row.center=TRUE,
	row.scale=TRUE, 
	low.threshold=-5, 
	high.threshold=5, 
	keep.jmat=FALSE,
	do.par = FALSE,
	num.cores = 1
){
	# step 0) check all the input parameters
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("'object' is not a snap object")
	}
			
	# 3. check if mat is binary;
	mat = match.arg(mat);
	if(mat == "bmat"){
		x = object@bmat;		
		if(nrow(x) == 0){
			stop("'bmat' is empty")		
		}
		if(max(object@bmat) > 1 || min(object@bmat) < 0){
			stop("'bmat' is not a binary matrix, run 'makeBinary' first")	
		}
	}else if(mat == "pmat"){
		x = object@pmat;
		if(nrow(x) == 0){
			stop("'pmat' is empty")		
		}
		if(max(object@pmat) > 1 || min(object@pmat) < 0){
			stop("'pmat' is not a binary matrix, run 'makeBinary' first")	
		}				
	}else if(mat == "gmat"){
		x = object@gmat;
		if(nrow(x) == 0){
			stop("'gmat' is empty")		
		}
		if(max(object@gmat) > 1 || min(object@gmat) < 0){
			stop("'gmat' is not a binary matrix, run 'makeBinary' first")	
		}				
	}else{
		stop("mat does not exist in object")
	}
	
	# 4. check if empty rows exist in bmat;
	if(any(Matrix::rowSums(x) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")	
	}
		
    # input checking for parallel options
    if (do.par) {
      if (num.cores == 1) {
        num.cores = 1
      } else if (num.cores > detectCores()) {
        num.cores <- detectCores() - 1
        warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
      }
    } else if (num.cores != 1) {
      num.cores <- 1
      warning("For parallel processing, please set do.par to TRUE.")
    }
	
	norm.method = match.arg(norm.method);
	# max.var can not exceed number of cells
	max.var = min(max.var, nrow(x));
		
	# step 1) randomly select a subset of cells as reference 
	obj.ref = object[sort(sample(seq(nrow(x)), max.var)),]
	
	# step 2) slice the orginal object into list
	id = seq(nrow(x));
	id.ls = split(id, ceiling(seq_along(id)/ncell.chunk));
	
	# always combine the last one with 2nd last one
	# just to prevent from having a chunk with less 
	# than 2 cells
	if(length(id.ls) > 1){
		id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
		# remove the last item of the list
		id.ls = id.ls[-length(id.ls)];
	}	
	
	obj.ls <- lapply(id.ls, function(x){
		object[x,]
	})
	
	message("Epoch: calculating jaccard index ...");
	# step 3) calculate jaccard index for each small pieces	
	obj.ls <- parallel::mclapply(obj.ls, function(x){
		calJaccard.snap(x, obj.ref, mat)
	}, mc.cores=num.cores);
	
	# step 4) normalize for each object chunk
	message("Epoch: normalizing jaccard index ...");
	obj.ls <- parallel::mclapply(obj.ls, function(x){
		normJaccard.snap(x, obj.ref, mat, method=norm.method, k)
	}, mc.cores=num.cores);
	
	message("Epoch: scaling normalized jaccard index ...");
	nmat <- do.call(rbind, lapply(obj.ls, function(x){x@nmat}));
	if(row.center || row.scale){
		nmat = t(scale(t(nmat), center=row.center, scale=row.scale));
	}
	
	# capping
	nmat[nmat >= high.threshold] = high.threshold;
	nmat[nmat <= low.threshold]  = low.threshold;
	
	# step 5) combine
	object@nmat = nmat;
	if(keep.jmat){
		object@jmat <- do.call(rbind, lapply(obj.ls, function(x){x@jmat}));
	}
	return(object);
}

calJaccard.snap <- function(obj1, obj2, mat=c("bmat", "pmat", "gmat")){
	mat = match.arg(mat);
	if(mat == "bmat"){
		X_i = obj1@bmat;
		X_j = obj2@bmat;
		A = Matrix::tcrossprod(X_i, X_j);
		bi = Matrix::rowSums(X_i);
		bj = Matrix::rowSums(X_j);
		obj1@jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));			
	}else if(mat=="pmat"){
		X_i = obj1@pmat;
		X_j = obj2@pmat;
		A = Matrix::tcrossprod(X_i, X_j);
		bi = Matrix::rowSums(X_i);
		bj = Matrix::rowSums(X_j);
		obj1@jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));					
	}else if(mat=="gmat"){
		X_i = obj1@gmat;
		X_j = obj2@gmat;
		A = Matrix::tcrossprod(X_i, X_j);
		bi = Matrix::rowSums(X_i);
		bj = Matrix::rowSums(X_j);
		obj1@jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));					
	}
	return(obj1);
}

normJaccard.snap <- function(obj1, obj2, mat=c("bmat", "pmat", "gmat"), method=c("normOVE", "normOVN"), k=15){
	mat = match.arg(mat);
	method = match.arg(method);
	if(mat == "bmat"){
		b1 = Matrix::rowMeans(obj1@bmat);
		b2 = Matrix::rowMeans(obj2@bmat);
		jmat = obj1@jmat;
		jmat[jmat == 1] = mean(jmat);
		x = jmat;	
		emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));			
		if(method == "normOVE"){
			data = data.frame(x=c(emat), y=c(jmat))	
			model = lm(y ~ x, data)
			nmat = matrix(model$residuals, nrow(emat), ncol(emat));
		}else{
			nmat = jmat - emat;
		}
		obj1@nmat = nmat;				
	}else if(mat=="pmat"){
		b1 = Matrix::rowMeans(obj1@pmat);
		b2 = Matrix::rowMeans(obj2@pmat);
		jmat = obj1@jmat;
		jmat[jmat == 1] = mean(jmat);	
		x = jmat;
		emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));			
		if(method == "normOVE"){
			data = data.frame(x=c(emat), y=c(jmat))	
			model = lm(y ~ x, data)
			nmat = matrix(model$residuals, nrow(emat), ncol(emat));
		}else{
			nmat = jmat - emat;
		}
		obj1@nmat = nmat;				
	}else if(mat=="gmat"){
		b1 = Matrix::rowMeans(obj1@gmat);
		b2 = Matrix::rowMeans(obj2@gmat);
		jmat = obj1@jmat;
		jmat[jmat == 1] = mean(jmat);	
		x = jmat;
		emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));			
		if(method == "normOVE"){
			data = data.frame(x=c(emat), y=c(jmat))	
			model = lm(y ~ x, data)
			nmat = matrix(model$residuals, nrow(emat), ncol(emat));
		}else{
			nmat = jmat - emat;
		}
		obj1@nmat = nmat;				
	}
	
	return(obj1);
}

# estimate the expected jaccard index using OVN
.normOVN <- function(o, p1, p2, k){
	# sort the jaccard index matrix based on the coverage
	ind1 = order(p1);
	ind2 = order(p2);
	o_srt = as.matrix(o[ind1, ind2, drop=FALSE]);
	# calculate expected jaccard index
	mask_mat <- matrix(1, k, k);
	exp = focal(raster(as.matrix(o_srt)), mask_mat, mean, na.rm=TRUE, pad = T);
	ee = raster::as.matrix(exp)[order(ind1),order(ind2),drop=FALSE];
	return(ee)
}

# estimate the expected jaccard index using OVE
.normOVE <- function(o, p1, p2, k){
    pp = tcrossprod(p1, p2);
	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
	ee = pp/(ss - pp)
	return(ee)	
}

#calJaccard.default <- function(object, cell.pcore=1000, ncore=1, max.var=5000, seed.use=10){
#	# 1. check if object is a snap;
#	if(class(object) != "snap"){
#		stop("'object' is not a snap object")
#	}
#
#	# 2. check of bmat exists;
#	if(nrow(object@bmat) == 0){
#		stop("'bmat' is empty")		
#	}
#
#	# 3. check if bmat is really binary;
#	# this is not strict!
#	if(max(object@bmat) > 1 || min(object@bmat) < 0){
#		stop("'bmat' is not a binary matrix, run 'makeBinary' first")	
#	}
#
#	# 4. check if empty rows exist in bmat;
#	if(any(Matrix::rowSums(object@bmat) == 0)){
#		stop("'bmat' contains empty rows, remove it first")	
#	}
#	
#	max.var = min(max.var, nrow(object@bmat));
#	
#	# number of cells
#	ncell = nrow(object@bmat);
#	object@idx = sort(sample(1:ncell, max.var));	
#			
#	# calculate the jaccard index
#	if(ncell <= cell.pcore){
#		X = object@bmat;		
#		X_i = X;
#		X_j = X[object@idx,];
#	    A = Matrix::tcrossprod(X_i, X_j);
#		bi = Matrix::rowSums(X_i);
#		bj = Matrix::rowSums(X_j);
#		jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));	
#	}else{
#		jmat = matrix(0, ncell, ncell);
#		# first cut the entire binary matrix into small pieces
#		id = 1:ncell;
#		id.ls = split(id, ceiling(seq_along(id)/cell.pcore));
#	
#		# always combine the last one with 2nd last one
#		# just to prevent from having a chunk with less 
#		# than 2 cells
#		if(length(id.ls) > 1){
#			id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
#			# remove the last item of the list
#			id.ls = id.ls[-length(id.ls)];
#		}
#	
#		# conquar each small problem sperately and combine them together
#		X = object@bmat;
#		jmat.chunk <- mclapply(id.ls, function(x){
#			id_i = x;
#			id_j = object@idx;
#			X_i = X[id_i,];
#			X_j = X[id_j,];
#		    A = Matrix::tcrossprod(X_i, X_j);
#			bi = Matrix::rowSums(X_i);
#			bj = Matrix::rowSums(X_j);
#			J_ij = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));			
#		}, mc.cores=ncore);		
#		jmat = do.call(rbind, jmat.chunk);
#	}		
#	
#	if(max(jmat) != 1 || min(jmat) < 0){
#		stop("@jmat error detected, jmat has values greater than 1 or smaller than 0, make sure @bmat is binary first")
#	}
#	# check if 
#	if(nrow(jmat) != nrow(object@bmat)){
#		stop("@jmat error detected, jmat is incomplete, try to use fewer cores!")		
#	}
#	object@jmat = jmat;
#	return(object);
#}
#
#calJaccard.default <- function(object, cell.pcore=1000, ncore=1, max.var=5000, seed.use=10){
#	# 1. check if object is a snap;
#	if(class(object) != "snap"){
#		stop("Not a snap object")
#	}
#		
#	# 2. check of bmat exists;
#	if(nrow(object@bmat) == 0){
#		stop("@bmat is empty")		
#	}
#
#	# 3. check if bmat is really binary;
#	# this is not strict!
#	if(max(object@bmat) > 1 || min(object@bmat) < 0){
#		stop("@bmat is not a binary matrix")	
#	}
#
#	# 4. check if empty rows exist in bmat;
#	if(any(Matrix::rowSums(object@bmat) == 0)){
#		stop("@bmat contains empty rows, remove it first")	
#	}
#
#	# 5. check of jmat already exists;
#	if(nrow(object@jmat) > 0){
#		message("@jmat already exists")	
#	}
#	
#	max.var = min(max.var, nrow(object@bmat));
#
#	# number of cells
#	ncell = nrow(object@bmat);
#	object@idx = sort(sample(1:ncell, max.var));	
#	
#	# calculate the jaccard index
#	if(ncell <= cell.pcore){
#		jmat = calJaccard.snap(object, object)
#	}else{
#		# first cut the entire binary matrix into small pieces
#		id = 1:ncell;
#		id.ls = split(id, ceiling(seq_along(id)/cell.pcore));
#		
#		# always combine the last one with 2nd last one
#		# just to prevent from having a chunk with less 
#		# than 2 cells
#		if(length(id.ls) > 1){
#			id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
#			# remove the last item of the list
#			id.ls = id.ls[-length(id.ls)];
#		}
#		
#		
#		# conquar each small problem sperately and combine them together
#		jmat = do.call(rbind, jmat.chunk);
#		}	
#}
#
#
#calJaccard.snap <- function(obj1, obj2){
#	jmat = calJaccard.bmat(obj1@bmat, obj2@bmat)
#	return(jmat);
#}
#


