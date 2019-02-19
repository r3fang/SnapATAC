#' Dimentionality Reduction for Visualization
#'
#' Wrapper for the implementation of Barnes-Hut t-Distributed
#' Stochastic Neighbor Embedding and FFT-accelerated Interpolation-based 
#' t-SNE (FIt-SNE) t-SNE, and UMAP
#'
#' @param object A snap object.
#' @param dims integer; Output dimensionality (default: 2)
#' @param pca_dims integer vector; dimentions of PCA used for tsne or umap
#' @param method character; A character variable indicates what t-sne method to use ("Rtsne", "fast_tsne", "umap").
#' @param fast_tsne_path character; A character variable indicates where the excutable fast_tsne located. Required when method="fast_tsne"
#' @param Y_init matrix; matrix; Initial locations of the objects. If NULL, random initialization will be used (default: NULL). Note that when using this, the initial stage with exaggerated perplexity values and a larger momentum term will be skipped.
#' @param seed.use numeric; random seed
#' @param tsne_perplexity numeric; perplexity for runing tsne
#' @param tsne_theta numeric; theta for tsne
#' @param max_tsne_iter numeric; max iteration for running tsne
#'
#' @return Returns a snap object with the tsne or umap
#'
#' @export

runViz<- function(object, ...) {
  UseMethod("runViz");
}

runViz.default <- function(
	object, 
	dims=2,
	pca_dims=NULL, 
	method=c("Rtsne", "umap", "fast_tsne"), 
	fast_tsne_path = NULL, 
	Y_init=NULL,
	seed.use=10, 
	tsne_perplexity=30, 
	tsne_theta = 0.5,
	max_tsne_iter=1000,
	...){
	# check input
	if(class(object) != "snap"){
		stop("'object' is not a snap object")
	}
	
	if(ncol(object@smat) == 0){
		stop("run 'runPCA' first")
	}
	# check PCA dimentions
	if(is.null(pca_dims)){
		pca_dims = seq(ncol(object@smat));
	}else{
		if(any((pca_dims %in% seq(ncol(object@smat))) == FALSE)){
			stop("'pca_dims' exceeds the PCA space, reset 'pca_dims' and run it again")
		}		
	}
	# check input parameters
	method = match.arg(method);
	# check if fi method exists or not
	if(method=="fast_tsne"){
		if(!file.exists(fast_tsne_path)){
			stop("'fast_tsne_path' fast tsne does not exist")
			set.seed(seed.use)
			object@tsne = fftRtsne(
				X = object@smat[,pca_dims], 
				dims=dims, 
				fast_tsne_path=fast_tsne_path, 
				perplexity=tsne_perplexity, 
				max_tsne_iter=max_tsne_iter, 
				rand_seed=seed.use);					
		}
	}else if(method=="Rtsne"){
		set.seed(seed.use)
		object@tsne = Rtsne(
			object@smat[,pca_dims], 
			dims=dims, 
			perplexity=tsne_perplexity, 
			max_iter = max_tsne_iter, 
			verbose = FALSE, 
			pca = FALSE, 
			is_distance = FALSE, 
			check_duplicates=FALSE)$Y;			
	}else{
		set.seed(seed.use)
		object@umap = umap(object@smat[,pca_dims])$layout;
	}
	return(object);
}
