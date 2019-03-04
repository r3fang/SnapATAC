#' Find Clusters Using Louvain Algorithm
#'
#' @param obj A snap object
#' @param k Neighbour size for creating K-nearest neighbour (knn) graph [k=15].
#' @param pca_dims Selected principal components to use for creating kNN graph.
#' @param method A character class that indicates what cluster method to use c("louvain", "jaccard_louvain") ["louvain"].
#' @param resolution A numeric class that indicates the resolution for louvain clustering [1].
#' @param data_path Directory path that stores temporary data [NULL].
#' @param result_path Directory path that stores clustering result [NULL].
#' @param path_to_snaptools Path to snaptools excutable binary file [NULL].
#' @return Returns a Snap obj with the cluster stored in obj@cluster
#'
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adjacency_matrix cluster_louvain
#' @importFrom utils file_test write.table read.table file.remove
#' 
#' @export

runCluster <- function(obj, k, pca_dims, method, resolution, data_path, result_path, path_to_snaptools) {
  UseMethod("runCluster", obj);
}

#' @export
runCluster.default <- function(
	obj, 
	k=15, 
	pca_dims=NULL, 
	resolution=1, 
	data_path=NULL, 
	method="louvain",
	result_path=NULL,
	path_to_snaptools=NULL
	){		
		
	if (is.null(data_path)) {
		data_path <- tempfile(pattern='community_data_', fileext='.dat')
	}
	
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='community_result_', fileext='.dat')
	}
	
	if(!is.null(path_to_snaptools)){
		path_to_snaptools <- normalizePath(path_to_snaptools);
		if (!file_test('-x', path_to_snaptools)) {
			warning(path_to_snaptools, " does not exist or is not executable; switch to igraph cluster_louvain function, but resolution is not supported")
			method = "R-igraph";			
		}else{
			method = "pyhton-louvain";			
		}		
	}else{
		method = "R-igraph";
		warning("Using igraph cluster_louvain function, resolution is not supported");
	}
	
	# 1. check if obj is a snap obj
	if(class(obj) != "snap"){
		stop("'obj' is not a snap obj");
	}

	# 2. check if smat exists
	if(nrow(obj@smat) == 0){
		stop("dimentionality reduction does not exist in obj, run runDimReduct first");
	}
	
	ncell = nrow(obj@smat);
	nvar = ncol(obj@smat);
	
	# 3. check pc.sel exceed the smat dimetions 
	if(is.null(pca_dims)){
		pca_dims=1:nvar;	
	}else{
		if(any(pca_dims > nvar) ){
			stop("'pca_dims' exceeds reduced dimentions variables number");
		}		
	}

	dx = nn2(obj@smat[,pca_dims], k = k+1)$nn.idx
	dx = dx[,2:(k+1)]
	edges <- matrix(unlist(sapply(1:nrow(dx),function(i) { rbind(rep(i,k),dx[i,])})),nrow=2);
	edges = as.data.frame(t(edges));	
	if(method == "R-igraph"){
		ajm = Matrix(0, ncell, ncell, sparse=TRUE);
		ajm[as.matrix(edges[,1:2])] = 1;
		g = graph_from_adjacency_matrix(ajm, mode="undirected", weighted=NULL);		
		rm(ajm);
		cl = cluster_louvain(g);
		obj@cluster = factor(cl$membership);	
	}else if(method == "python-louvain"){
		edges$weight = 1;
		write.table(edges, file = data_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
		flag = system2(command=path_to_snaptools, args=c("louvain", "--edge-file", data_path, "--output-file", result_path, "--resolution", resolution));		
		if (flag != 0) {
		   	stop("'runCluster' call failed");
		}	
		cluster = read.table(result_path);		
		file.remove(data_path)
		file.remove(result_path)
		obj@cluster = factor(cluster[order(cluster[,1]),2]);	
	}
	return(obj);
}

