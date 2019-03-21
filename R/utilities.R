
# Find centroid of each cluster
# @param x coordinates
# @param y cluster label
findCentrod <- function(x, y){
	x.ls = split(data.frame(x),y);
	centroid.ls = lapply(split(data.frame(x),y), colMeans)
	centroid.df = data.frame(do.call(rbind, centroid.ls))
	centroid.df$Y = names(centroid.ls);
	return(centroid.df);		
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


calJaccardSingle <- function(file.name, mat_j){
	mat_i = readRDS(file.name);
	calJaccard(mat_i, mat_j);
}

# Normlize jaccard index
# @param jmat Jaccard index matrix
# @param b1 Coverage probability for rows
# @param b2 Coverage probability for columns
# @param method Method used for normalization c("normOVE", "normOVN")
# @param k neighbour size for normOVN
normJaccard <- function(jmat, b1, b2, method, k=15){
	# estimate the expected jaccard index using OVN
	#' @importFrom raster focal raster
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
	
	jmat[jmat == 1] = mean(jmat);
	x = jmat;
	emat = eval(parse( text = paste(".", method, "(x, b1, b2, k)", sep="")));				
	if(method == "normOVE"){
		data = data.frame(x=c(emat), y=c(jmat))	
		model = stats::lm(y ~ x, data)
		nmat = matrix(model$residuals, nrow(emat), ncol(emat));
	}else if(method == "normOVN"){
		nmat = jmat - emat;
	}
	rm(jmat)		
	rm(emat)		
	return(nmat);
}

findNegCells <- function(obj, idx.pos, method=c("knn", "random", "other")){
	method=match.arg(method);
	mat.use = getGraph(obj@graph);
	ncell = nrow(obj);
	ncell.pos = length(idx.pos);
	ncell.neg = min(length(idx.pos), ncell - ncell.pos);
	
	if(method == "knn"){
		idx.neg = setdiff(seq(ncell), idx.pos);
		dx = order(Matrix::colSums(mat.use[idx.pos,idx.neg]), decreasing=TRUE)[1:ncell.neg];
		idx.neg = idx.neg[dx]
	}else if(method == "random"){
		idx.neg = setdiff(seq(ncell), idx.pos);
		nn.num = min(length(idx.pos), length(idx.neg));
		set.seed(10);
		idx.neg = sample(idx.neg, nn.num);	
	}else if(method == "other"){
		idx.neg = setdiff(seq(ncell), idx.pos);			
	}
	return(idx.neg);
}

pvalue2fdr <- function(p1, p2){
	fdr.ls <- lapply(as.list(seq(0.001, 1, by=0.001)), function(p_i){
		(length(which(p2 <= p_i))) / (length(which(p1 <= p_i)))
	})
	fdr.tab = data.frame(p=seq(0.001, 1, by=0.001), fdr=do.call(c, fdr.ls));	
	return(fdr.tab);
}


#' @importFrom parallel mclapply
splitBmat <- function(mat, id.ls, num.cores=1){
	fileList = lapply(id.ls, function(x){
		tempfile(fileext=".rds")
	})
	mclapply(as.list(seq(id.ls)), function(i){
		saveRDS(mat[id.ls[[i]],], file=fileList[[i]]);
	}, mc.cores=num.cores);
	return(fileList);
}

