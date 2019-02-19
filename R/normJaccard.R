#' Normalize Jaccard Index Matrix
#'
#' This function takes a Snap object as input with jmat attributes and normalize 
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
#' @param object A snap object
#' @param method A character class that indicates the normalization method to be used. This must be one of "OVN", "OVE"
#' @param k A numeric class - number of neibouring cells to use for OVN (used only if method="OVN") (default k = 15) 
#' @param max.cell A numeric class that indicates the max number of cells to calculate per core
#' @param cores A numeric class that indicates the number of cores to use for calculation
#' @param outlier.filter A numeric class that indicates the top outlier.filter of values as outliers to be capped (0.01)
#' @param row.norm a logical value indicating whether rows of the normalized jaccard index matrix should be centered to 0.
#'
#' @export
#'
#normJaccard <- function(object, ...) {
#  UseMethod("normJaccard");
#}

#normJaccard.default <- function(object, k=15, outlier.filter=1e-3, eps=1e-6, cell.pcore=1000, ncore=5, method=c("normOVE", "normOVN"), row.center=TRUE, row.scale=TRUE){	
#	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
#	if (!is.wholenumber(k) || k<=0) { stop("Incorrect k (OVN).")}
#	if (!is.numeric(outlier.filter) || outlier.filter<=0 || outlier.filter > 1) { stop("Incorrect outlier.filter.")}
#	
#	method <- match.arg(method);	
#	
#	# 1. check if object is a snap;
#	if(class(object) != "snap"){
#		stop("Not a snap object")
#	}
#		
#	# 2. check of jmat exists;
#	if(nrow(object@jmat) == 0){
#		stop("@jmat is empty")		
#	}	
#
#	# 3. check of emat exists;
#	if(nrow(object@nmat) != 0){
#		message("Warning: @nmat already exists")		
#	}	
#	
#	# 3. make sure k is odd number
#	if((k %% 2) == 0) {
#		k = k + 1
#	}
#	
#	# 4. check idx
#	if(length(object@idx) == 0){
#		stop("@idx is empty");		
#	}else{
#		if(any(!object@idx %in% 1:nrow(object))){			
#			stop("@idx exceeds cell number");		
#		}
#	}
#	
#	row.means = Matrix::rowMeans(object@bmat);
#	ncell = length(row.means);
#	
#	# number of cells
#	if(ncell <= cell.pcore){
#		x = object@jmat;
#		x[x == 1] = mean(x);		
#		b1 = row.means;
#		b2 = row.means[object@idx];
#		emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));	
#		x[x == 1] = mean(x)
#		data = data.frame(x=c(emat), y=c(x))
#		model = lm(y ~ x, data)
#		nmat = matrix(model$residuals, nrow(emat), ncol(emat));
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
#		jmat = object@jmat
#		# conquar each small problem sperately and combine them together
#		nmat.chunk <- mclapply(id.ls, function(id_i){	
#			id_j = object@idx;
#			x = jmat[id_i,];
#			x[x == 1] = mean(x);		
#			b1 = row.means[id_i];
#			b2 = row.means[id_j];
#			emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));	
#			x[x == 1] = mean(x)
#			data = data.frame(x=c(emat), y=c(x))
#			model = lm(y ~ x, data)
#			nmat = matrix(model$residuals, nrow(emat), ncol(emat));
#		}, mc.cores=ncore);		
#		nmat = do.call(rbind, nmat.chunk);
#	}
#	
#	# scale emat to the same level with jaccard index
#	#y = object@jmat
#	#y[y == 1] = mean(y)
#	#data = data.frame(x=c(emat), y=c(y))
#	#model = lm(y ~ x, data)
#	#nmat = matrix(model$residuals, nrow(emat), ncol(emat));
#
#	# if row.norm is true
#	if(row.center || row.scale){
#		nmat = t(scale(t(nmat), center=row.center, scale=row.scale));
#	}
#	
#	# check if nmat is complete 
#	if(nrow(nmat) != nrow(object@bmat)){
#		stop("@nmat error detected, nmat is incomplete, try to use fewer cores!")		
#	}
#	object@nmat = nmat;
#	return(object);
#}
#
## estimate the expected jaccard index using OVN
#.normOVN <- function(o, p1, p2, k){
#	# sort the jaccard index matrix based on the coverage
#	ind1 = order(p1);
#	ind2 = order(p2);
#	
#	o_srt <- as.matrix(o[ind1, ind2, drop=FALSE]);
#    
#	# calculate expected jaccard index
#	mask_mat <- matrix(1, k, k);
#	exp = focal(raster(as.matrix(o_srt)), mask_mat, mean, na.rm=TRUE, pad = T);
#	ee = raster::as.matrix(exp)[order(ind1),order(ind2),drop=FALSE];
#	return(ee)
#}
#
## estimate the expected jaccard index using OVE
#.normOVE <- function(o, p1, p2, k){
#    pp = tcrossprod(p1, p2);
#	ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
#	ee = pp/(ss - pp)
#	return(ee)	
#}
