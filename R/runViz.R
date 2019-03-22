#' Dimentionality Reduction for Visualization
#'
#' Wrapper for the implementation of Barnes-Hut t-Distributed
#' Stochastic Neighbor Embedding and FFT-accelerated Interpolation-based 
#' t-SNE (FIt-SNE) t-SNE, and UMAP
#'
#' @param obj A snap object.
#' @param dims integer; Output dimensionality (default: 2)
#' @param pca.dims integer vector; principal components used for running tsne or umap
#' @param weight.by.sd Weight the cell embeddings by the sd of each PC
#' @param method character; A character variable indicates what t-sne method to use ("Rtsne", "fast_tsne", "umap").
#' @param fast_tsne_path character; A character variable indicates the path of the excutable FIt-sen located. Required only when method="fast_tsne"
#' @param Y.init matrix; Initial locations of the objs. If NULL, random initialization will be used (default: NULL). Note that when using this, the initial stage with exaggerated perplexity values and a larger momentum term will be skipped.
#' @param seed.use number; Random seeds to use.
#' @param num.cores number; Number of CPU used for computing.
#' @param tmp.folder Directory to store temporary files.
#' @param ... Arguments passed to Rtsne, umap or FIt-tsne
#' 
#' @return Returns a snap obj with the visulization
#' @importFrom Rtsne Rtsne
#' @importFrom utils installed.packages
#' @export

runViz<- function(obj, dims, pca.dims, weight.by.sd, method, fast_tsne_path, Y.init, seed.use, num.cores, tmp.folder, ...) {
  UseMethod("runViz", obj);
}

#' @export
runViz.default <- function(
	obj, 
	dims=2,
	pca.dims=NULL, 
	weight.by.sd=FALSE,
	method=c("Rtsne", "umap", "fast_tsne"), 
	fast_tsne_path = NULL, 
	Y.init=NULL,
	seed.use=10,
	num.cores=1,
	tmp.folder,
	...){
		
	# check input
	if(!is(obj, "snap")){
		stop("obj is not a snap object")
	}
	
	if(nrow(obj@smat@dmat) == 0){
		stop("PCA is empty, runDimReduct first")
	}

	if(nrow(obj@smat@dmat) != nrow(obj)){
		stop("PCA has different length with obj, data has been subsetted by mistake")
	}

	
	# check PCA dimentions
	ncell = nrow(obj);
	nvar = ncol(obj@smat@dmat);	
	if(is.null(pca.dims)){
		pca.dims=1:nvar;	
	}else{
		if(any(pca.dims > nvar) ){
			stop("'pca.dims' exceeds reduced dimentions variables number");
		}		
	}
	
	if(!is.logical(weight.by.sd)){
		stop("weight.by.sd must be a logical variable")
	}
	
	if(weight.by.sd){
		data.use = obj@smat@dmat[,pca.dims] %*% diag(sqrt(obj@smat@sdev[pca.dims])) ;
	}else{
		data.use = obj@smat@dmat[,pca.dims];
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
	
	# check input parameters
	method = match.arg(method);
	# check if fi method exists or not
	if(method=="fast_tsne"){
		if(is.null(fast_tsne_path)){
			stop("fast_tsne_path is missing");
		}

		if(!file.exists(fast_tsne_path)){
			stop("'fast_tsne_path' fast tsne does not exist")
		}		

		fast_tsne_path <- normalizePath(fast_tsne_path);
		if (!file_test('-x', fast_tsne_path)) {
			stop(fast_tsne_path, " is not executable; check your fast_tsne_path parameter")
		}				
		obj@tsne = fftRtsne(
				X = data.use, 
				dims=dims, 
				fast_tsne_path=fast_tsne_path, 
				rand_seed=seed.use,
				nthreads=num.cores,
				initialization=Y.init,
				...
			);	
		colnames(obj@tsne) = c("tsne-1", "tsne-2")				
	}else if(method=="Rtsne"){
		set.seed(seed.use);
		obj@tsne = Rtsne(
			data.use, 
			dims=dims, 
			verbose = FALSE, 
			pca = FALSE, 
			is_distance = FALSE, 
			check_duplicates=FALSE,
			num_threads=num.cores,
			Y_init=Y.init,
			rand_seed=seed.use,
			...
			)$Y;
		colnames(obj@tsne) = c("tsne-1", "tsne-2")
	}else{
		if (requireNamespace("umap", quietly = TRUE)) {
				set.seed(seed.use);
				obj@umap = umap::umap(data.use)$layout;
				colnames(obj@umap) = c("umap-1", "umap-2")
		  } else {
		      stop("Please install umap - learn more at https://cran.r-project.org/web/packages/umap/index.html")
		  }
		  
	}
	return(obj);
}


fftRtsne <- function(X, 
		     dims=2, perplexity=30, theta=0.5,
		     check_duplicates=TRUE,
		     max_iter=1000,
		     fft_not_bh = TRUE,
		     ann_not_vptree = TRUE,
		     stop_early_exag_iter=250,
		     exaggeration_factor=12.0, no_momentum_during_exag=FALSE,
		     start_late_exag_iter=-1.0,late_exag_coeff=1.0,
             mom_switch_iter=250, momentum=.5, final_momentum=.8, learning_rate=200,
		     n_trees=50, search_k = -1,rand_seed=-1,
		     nterms=3, intervals_per_integer=1, min_num_intervals=50, 
		     K=-1, sigma=-30, initialization=NULL,
		     data_path=NULL, result_path=NULL,
		     load_affinities=NULL,
		     fast_tsne_path=NULL, nthreads=0, perplexity_list = NULL, 
             get_costs = FALSE, df = 1.0,
			 tmp.folder,
			 ... ) {
        version_number = '1.1.0'

	if (is.null(fast_tsne_path)) {
		stop("fast_tsne_path is NULL")
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
	
	if (is.null(data_path)) {
		data_path <- tempfile(pattern='fftRtsne_data_', tmpdir = tmp.folder, fileext='.dat')
	}
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='fftRtsne_result_', tmpdir = tmp.folder, fileext='.dat')
	}
	if (is.null(fast_tsne_path)) {
		fast_tsne_path <- system2('which', 'fast_tsne', stdout=TRUE)
	}

	fast_tsne_path <- normalizePath(fast_tsne_path)
	if (!file_test('-x', fast_tsne_path)) {
		stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
	}

	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	if (!is.numeric(theta) || (theta<0.0) || (theta>1.0) ) { stop("Incorrect theta.")}
	if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
	if (!is.matrix(X)) { stop("Input X is not a matrix")}
	if (!(max_iter>0)) { stop("Incorrect number of iterations.")}
	if (!is.wholenumber(stop_early_exag_iter) || stop_early_exag_iter<0) { stop("stop_early_exag_iter should be a positive integer")}
	if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
	if (!is.numeric(df)) { stop("df should be numeric")}
	if (!is.wholenumber(dims) || dims<=0) { stop("Incorrect dimensionality.")}
	if (search_k == -1) {
       if (perplexity>0) {
          search_k = n_trees*perplexity*3
       } else if (perplexity==0) {
          search_k = n_trees*max(perplexity_list)*3
       } else { 
          search_k = n_trees*K
       }
    }

	if (fft_not_bh){
	  nbody_algo = 2;
	}else{
	  nbody_algo = 1;
	}

	if (is.null(load_affinities)) {
		load_affinities = 0;
	} else {
		if (load_affinities == 'load') {
			load_affinities = 1;
		} else if (load_affinities == 'save') {
			load_affinities = 2;
		} else {
			load_affinities = 0;
		}
	}
	
	if (ann_not_vptree){
	  knn_algo = 1;
	}else{
	  knn_algo = 2;
	}
	tX = c(t(X))

	f <- file(data_path, "wb")
	n = nrow(X);
	D = ncol(X);
	writeBin(as.integer(n), f,size= 4)
	writeBin( as.integer(D),f,size= 4)
	writeBin( as.numeric(theta), f,size= 8) #theta
	writeBin( as.numeric(perplexity), f,size= 8) #theta

    if (perplexity == 0) {
    	writeBin( as.integer(length(perplexity_list)), f, size=4)
	    writeBin( perplexity_list, f) 
    }

	writeBin( as.integer(dims), f,size=4) #theta
	writeBin( as.integer(max_iter),f,size=4)
	writeBin( as.integer(stop_early_exag_iter),f,size=4)
	writeBin( as.integer(mom_switch_iter),f,size=4)
	writeBin( as.numeric(momentum),f,size=8)
	writeBin( as.numeric(final_momentum),f,size=8)
	writeBin( as.numeric(learning_rate),f,size=8)
	writeBin( as.integer(K),f,size=4) #K
	writeBin( as.numeric(sigma), f,size=8) #sigma
	writeBin( as.integer(nbody_algo), f,size=4)  #not barnes hut
	writeBin( as.integer(knn_algo), f,size=4) 
	writeBin( as.numeric(exaggeration_factor), f,size=8) #compexag
	writeBin( as.integer(no_momentum_during_exag), f,size=4) 
	writeBin( as.integer(n_trees), f,size=4) 
	writeBin( as.integer(search_k), f,size=4) 
	writeBin( as.integer(start_late_exag_iter), f,size=4) 
	writeBin( as.numeric(late_exag_coeff), f,size=8) 
	
	writeBin( as.integer(nterms), f,size=4) 
	writeBin( as.numeric(intervals_per_integer), f,size=8) 
	writeBin( as.integer(min_num_intervals), f,size=4) 
	tX = c(t(X))
	writeBin( tX, f) 
	writeBin( as.integer(rand_seed), f,size=4) 
        writeBin(as.numeric(df), f, size=8)
	writeBin( as.integer(load_affinities), f,size=4) 
	if (! is.null(initialization)){ writeBin( c(t(initialization)), f) }		
        print(df)
	close(f) 

	flag= system2(command=fast_tsne_path, args=c(version_number,data_path, result_path, nthreads));
	if (flag != 0) {
		stop('tsne call failed');
	}
	f <- file(result_path, "rb")
	n <- readBin(f, integer(), n=1, size=4);
	d <- readBin(f, integer(), n=1, size=4);
	Y <- readBin(f, numeric(), n=n*d);
        Y <- t(matrix(Y, nrow=d));
        if (get_costs ) {
            tmp <- readBin(f, integer(), n=1, size=4);
            costs <- readBin(f, numeric(), n=max_iter,size=8);
            Yout <- list( Y=Y, costs=costs);
        }else{
            Yout <- Y;
        }
        close(f)
        file.remove(data_path)
        file.remove(result_path)
        return(Yout)
}




