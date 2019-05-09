#' Latent Dirichlet Allocation
#'
#' This function takes a Snap obj as input with bmat/pmat slot and run 
#' Latent Dirichlet Allocation (LDA).  
#' 
#' LDA is a generative statistical model that allows sets of observations to be explained by unobserved groups 
#' that explain why some parts of the data are similar. was first applied to analyze single cell ATAC-seq by 
#' Cis-Topic (González-Blas, Nature Methods, 2019). LDA iteratively optimize two probability distributions: 
#' (1) the probability of a region belonging to a topic (region–topic distribution) and 
#' (2) the contribution of a topic within a cell (topic–cell distribution).
#' 
#' Multiple LDA models will be trained and the optimal one will be selected according to the likelohood.
#' 
#' @param obj A snap obj
#' @param input.mat Input matrix to be used for LSA c("bmat", "pmat").
#' @param topic An integer number of topics to return c(10, 20, 30).
#' @param method Method used to normalize the cell-by-topic score c("Z-score", "Probability").
#' @param num.cores Number of cores used for computing [1].
#' @param min.cell Min cell coverage. Features with coverage less than min.cell will be filtered [10].
#' @param iterations The number of sweeps of Gibbs sampling over the entire corpus to make.
#' @param burnin A scalar integer indicating the number of Gibbs sweeps to
#'	consider as burn-in (i.e., throw away) for
#'	'lda.collapsed.gibbs.sampler' and
#'	'mmsb.collapsed.gibbs.sampler'.  If this parameter is
#'	non-NULL, it will also have the side-effect of enabling the
#'	document_expects field of the return value (see below for
#'	details).  Note that burnin iterations do NOT count towards
#'	num.iterations.
#' @param alpha The scalar value of the Dirichlet hyperparameter for topic proportions.
#' @param beta The scalar value of the Dirichlet hyperparamater for topic multinomials.
#' @param alphaByTopic scale alpha by topic number
#' @param seed.use A numeric class that indicates random seeding number [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runLDA(
#'	obj=demo.sp, 
#'	input.mat="bmat", 
#'	topic=c(10, 20, 30), 
#'	method="Z-score", 
#'	min.cell=0,
#'	num.cores=3
#'	);
#' 
#' @import Matrix
#' @importFrom parallel makeCluster clusterEvalQ
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom plyr laply 
#' @export

runLDA <- function(
	obj,
	input.mat = c("bmat", "pmat"),
	topic=c(10, 20, 30),
	method=c("Z-score", "Probability"),
	num.cores=1,
	min.cell=10,
	seed.use=10,
	iterations = 500,
	burnin = 250,
	alpha = 50,
	alphaByTopic = TRUE,
	beta=0.1
){
	message("Epoch: checking the inputs ...");
	# 1. check if lda has been installed
	if (!requireNamespace("lda", quietly = TRUE)) {
      stop("Please install R package lda first - install.packages(lda)")
	}
	
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
		data.use = t(obj@bmat);
	}else if(input.mat == "pmat"){
		data.use = t(obj@pmat);
	}else{
		stop("input.mat does not exist in obj")
	}
	
	if (burnin >= iterations){
		stop('The number of iterations must be higher than the burnin!')
	} 
  	
	if((x = max(data.use@x)) > 1L){
		stop("input matrix is not a binary matrix, run 'makeBinary' first");
	}
	
	ncell = ncol(data.use);
	if(!is(min.cell, "numeric")){
		stop("min.cell is not a numeric");
	}else{
		if(min.cell > ncell){
			stop("min.cell exceeds number of cells");
		}
	}
	
	message("Epoch: filtering low-coverage features ...");
	idx = which(Matrix::rowSums(data.use) > min.cell);
	if((x = length(idx)) == 0L){
		stop("No features has coverage than min.cell");
	}else{
		data.use = data.use[idx,];
	}
	
	cell.cov = Matrix::colSums(data.use);
	if(any(cell.cov == 0)){
		stop("empty cells detected after removing low-coverage features, remove empty cells first.")
	}
	
 	# Take binary count matrix
 	cellnames <- obj@barcode;
	regionnames <- paste("region", seq(nrow(data.use)), sep="."); 
	
	x =  Matrix::colSums(data.use);
	if(any((x==0))){
		stop("empty rows found in the input matrix");
	}
	
	message("Epoch: preparing the data ...");
    # Prepare data
    cellList <- split(
		as.integer(data.use@i), 
		rep(seq_along(data.use@p+1), 
		times=diff(c(data.use@p+1, length(data.use@i) + 1)))
	)
    rm(data.use);

    cellList <- lapply(cellList, function(x) {
		x <- rbind(x, rep(as.integer(1), length(x)))
	}) 

    names(cellList) <- cellnames
    cellList <- lapply(cellList, function(x) {colnames(x) <- regionnames[x[1,]+1];x})
    regionList <- regionnames
  	
	message("Epoch: running LDA models ...");
  	if (length(topic) > 1){
    	if (length(topic) < num.cores){
      	  print(paste('The number of cores (', num.cores, ') is higher than the number of models (', length(topic),').', sep=''))
    }
    	if (num.cores > 1){
    	  cl <- makeCluster(num.cores, type = "SOCK");
    	  registerDoSNOW(cl);
    	  clusterEvalQ(cl, library(lda));
    	  clusterExport(cl, c("cellList", "regionList", "topic", "iterations", "burnin", "alpha", "beta"), envir=environment());
    	  opts <- list(preschedule=TRUE)
    	  clusterSetRNGStream(cl, seed.use)
    	  if (alphaByTopic==TRUE){
    	    models <- suppressWarnings(
				plyr::llply(
					.data=topic, 
					.fun=function(t) 
						lda.collapsed.gibbs.sampler(
							documents=cellList, 
							K=t, 
							vocab=regionList, 
							num.iterations=iterations, 
							alpha=alpha/t, eta=beta, 
							compute.log.likelihood = TRUE, 
							burnin=burnin
							)[-1], 
						.parallel = TRUE, 
						.paropts = list(.options.snow=opts), 
						.inform = FALSE
					)
				)
    	  
		  }else{
    	    models <- suppressWarnings(
				plyr::llply(
					.data=topic, 
					.fun=function(t) 
						lda.collapsed.gibbs.sampler(
							cellList, 
							t, 
							regionList, 
							num.iterations=iterations, 
							alpha=alpha, eta=beta, 
							compute.log.likelihood = TRUE, 
							burnin=burnin
							)[-1], 
						.parallel = TRUE, 
						.paropts = list(.options.snow=opts), 
						.inform=FALSE
					)
				)
			
    	  }
    	  stopCluster(cl)
	  }else{
		  if (alphaByTopic==TRUE){
    	    models <- suppressWarnings(
				plyr::llply(
					.data=topic, 
					.fun=function(t) 
						lda.collapsed.gibbs.sampler(
							cellList, 
							t, 
							regionList, 
							num.iterations=iterations, 
							alpha=alpha/t, eta=beta, 
							compute.log.likelihood = TRUE, 
							burnin=burnin)[-1], 
						.progress = progress_text(char = ".")
			  		  )
		 		)
    	  }
    	  else{
      	    models <- suppressWarnings(
  				plyr::llply(
  					.data=topic, 
  					.fun=function(t) 
  						lda.collapsed.gibbs.sampler(
  							cellList, 
  							t, 
  							regionList, 
  							num.iterations=iterations, 
  							alpha=alpha, eta=beta, 
  							compute.log.likelihood = TRUE, 
  							burnin=burnin)[-1], 
  						.progress = progress_text(char = ".")
  					)
				)
    	  }
    }
	
	}else{
		set.seed(seed.use)
		if (alphaByTopic==TRUE){
			models <- plyr::llply(
				lda.collapsed.gibbs.sampler(
					cellList, topic, regionList, 
					num.iterations=iterations, 
					alpha=alpha/topic, eta=beta, 
					compute.log.likelihood = TRUE, 
					burnin=burnin)[-1], 
				.progress = progress_text(char = ".")
			)
		}else{
			models <- plyr::llply(
				lda.collapsed.gibbs.sampler(
					cellList, topic, regionList, 
					num.iterations=iterations, 
					alpha=alpha, eta=beta, 
					compute.log.likelihood = TRUE, 
					burnin=burnin)[-1], 
				.progress = progress_text(char = ".")
			)
		}
    }
	
	message("Epoch: selecting LDA model ...");
  	names(models) <- laply(1:length(models), function(x) sapply(models[x], function(y) nrow(y$topic_sums)));
  	loglikelihood <- sapply(seq_along(models), FUN=function(i) models[[i]]$log.likelihood[2,ncol(models[[i]]$log.likelihood)]);
	topics <-  sapply(seq_along(models), FUN=function(i) nrow(models[[i]]$topics));
	object.log.lik <- data.frame(topics=topics, LL=loglikelihood);
	object.log.lik <- object.log.lik[order(object.log.lik$topics),];	
	selected.model <- models[[which(object.log.lik$LL == max(object.log.lik$LL))]];
	
    if (method == 'Z-score'){
       modelMat <- scale(selected.model$document_expects, center=TRUE, scale=TRUE)
     }else if (method == 'Probability'){
      alpha <- alpha/length(selected.model$topic_sums)
      modelMat <- apply(selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }

	obj@smat@dmat = t(modelMat);
	obj@smat@sdev = rep(1, nrow(modelMat));
	obj@smat@iter = iterations;
	obj@smat@method = "LDA";				
	obj@smat@imat = input.mat;				
	message("Epoch: Done");
	return(obj)
}

