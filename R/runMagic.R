#' Affinity Graph-based Smoothing
#'
#' This function takes a snap obj as input with bmat/pmat/gmat slot and run 
#' magic to smooth the signal. 
#'
#' @param obj A snap obj
#' @param input.mat Input matrix c("bmat", "pmat", "gmat", "mmat").
#' @param step.size Number of steps for diffusion [3].
#'
#' @import Matrix
#' @export
runMagic <- function(
	obj,
	input.mat=c("bmat", "pmat", "gmat", "mmat"),
	step.size=3
){
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
		data.use = obj@mmat;	
		peak.use = GenomicRanges::GRanges();	
	}
		
	obj = obj[,idy.use, mat=input.mat];
	
	num.total.cells = nrow(data.use);
	if(num.landmarks >= num.total.cells){
		message("Epoch: computing jaccard similarity matrix ...");
		obj = runJaccard2(obj, obj);
		
		if(norm){
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
		obj@metaData$landmark = rep(1, num.total.cells);
	}else{
		message("Epoch: sampling landmarks ...");
		row.covs = log((Matrix::rowSums(data.use))+1,10);		
		row.covs.dens <- density(x = row.covs, bw = 'nrd', adjust = 1);
		sampling_prob <- 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps);
		set.seed(seed.use);
		idx.landmark <- sort(sample(x = seq(num.total.cells), size = num.landmarks, prob = sampling_prob));
		idx.other <- setdiff(seq(num.total.cells), idx.landmark);
		obj.landmark = obj[idx.landmark,]

		message("Epoch: computing jaccard similarity matrix ...");
		obj.landmark = runJaccard2(obj.landmark, obj.landmark);
		
		if(norm){
			message("Epoch: fitting regression model ...");
			obj.landmark = trainRegression(obj.landmark)
			message("Epoch: performing normalization ...");
			obj.landmark = normJaccard(obj.landmark, obj.landmark@regModel[1], obj.landmark@regModel[2], obj.landmark@regModel[3]);			
		}else{
			obj.landmark@jmat@nmat = obj.landmark@jmat@jmat;
		}

		# cap the max value
		nmat.cutoff = quantile(obj.landmark@jmat@nmat, nmat.outlier);
		obj.landmark@jmat@nmat[obj@jmat@nmat > nmat.cutoff] = nmat.cutoff;
		
		message("Epoch: computing eigen decomposition ...");
		obj.landmark = runEigDecomp(obj.landmark, num.eigs);

		idx.other.ls = split(idx.other, ceiling(seq_along(seq(idx.other))/ncell.chunk))
		eig.mat =list()
		# TODO - PARALLEL
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
		obj@metaData$landmark = rep(0, num.total.cells);
		obj@metaData$landmark[idx.landmark] = 1;
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
