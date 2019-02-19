#' Find Clusters Using Louvain Method
#'
#' This function takes a snap object with "smat" attributes and identify cluster using graph-based clustering method
#'
#' @param object A snap object.
#' @param k A numeric object indicates the size for KNN graphs [15].
#' @param pc.sel A numeric vector indicates which PCs are used for creating KNN graph [NULL].
#' @param func A character object indicates what louvain function to use "python" or "R" [python] .
#' @param method A character object indicates what cluster methods to use "louvain" or "jaccard.louvain" [louvain].
#' @param resolution A numeric vector between 0 and 1 indicates what resolution to use for louvain [0.5].
#' @param data_path A character variable indicates the path for storing graph data [NULL].
#' @param result_path A character variable indicates the path for storing clustering result [NULL].
#' @param community_path A character variable indicates the path to community_louvain program [bin/community_louvain].
#'
#' @return Returns a Snap object with the cluster stored in object@cluster
#'
#' @export
	
runCluster <- function(object, ...) {
  UseMethod("runCluster", object);
}

#' @export
runCluster.default <- function(
	object, 
	k=15, 
	centers=NULL,
	pca_dims=NULL, 
	method = c("jaccard_louvain", "louvain", "kmeans", "densityClust"), 
	resolution=1, 
	data_path=NULL, 
	result_path=NULL,
	path_to_louvain=NULL
){		
	if (is.null(data_path)) {
		data_path <- tempfile(pattern='community_data_', fileext='.dat')
	}
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='community_result_', fileext='.dat')
	}
	if(method %in% c("jaccard_louvain", "louvain")){
		if (is.null(path_to_louvain)) {
			path_to_louvain <- system2('which', 'community_louvain', stdout=TRUE)
		}
		path_to_louvain <- normalizePath(path_to_louvain)

		if (!file_test('-x', path_to_louvain)) {
			stop(path_to_louvain, " does not exist or is not executable; check your fast_tsne_path parameter")
		}		
	}
	
	method = match.arg(method);
	# 1. check if object is a snap object
	if(class(object) != "snap"){
		stop("'object' is not a snap object");
	}

	# 2. check if smat exists
	if(nrow(object@smat) == 0){
		stop("smat does not exist in object");
	}
	
	ncell = nrow(object@smat);
	nvar = ncol(object@smat);
	
	# 3. check pc.sel exceed the smat dimetions 
	if(is.null(pca_dims)){
		pca_dims=1:nvar;	
	}else{
		if(any(pca_dims > nvar) ){
			stop("'pca_dims' exceeds smat's variables number");
		}		
	}
	
	if(method == "jaccard_louvain"){
		dx = nn2(object@smat[,pca_dims], k = k+1)$nn.idx
		dx = dx[,2:(k+1)]
		edges <- matrix(unlist(sapply(1:nrow(dx),function(i) { rbind(rep(i,k),dx[i,])})),nrow=2);
		edges = as.data.frame(t(edges));	
		edges$weight = 1;
		ajm = Matrix(0, ncell, ncell, sparse=TRUE);
		ajm[as.matrix(edges[,1:2])] = edges[,3];
		g = graph_from_adjacency_matrix(ajm, mode="undirected", weighted=NULL);		
		g = graph_from_adjacency_matrix(similarity(g, method = "jaccard"), mode="undirected", weighted=TRUE)			
		adj=get.adjacency(g,attr='weight') 
		edges = as.data.frame(Matrix::which(adj > 0, arr.ind=TRUE))
		edges = edges[edges[,1] < edges[,2],]
		edges$value = adj[as.matrix(edges)]
		write.table(edges, file = data_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
		flag = system2(command=path_to_louvain, args=c("-i", data_path, "-o", result_path, "-r", resolution));
		if (flag != 0) {
		   	stop("'runCluster' call failed");
		}
	
		cluster = read.table(result_path);
		object@cluster = factor(cluster[order(cluster[,1]),2]);
		file.remove(data_path)
		file.remove(result_path)
	}else if(method == "louvain"){
		dx = nn2(object@smat[,pca_dims], k = k+1)$nn.idx
		dx = dx[,2:(k+1)]
		edges <- matrix(unlist(sapply(1:nrow(dx),function(i) { rbind(rep(i,k),dx[i,])})),nrow=2);
		edges = as.data.frame(t(edges));	
		edges$weight = 1;
		write.table(edges, file = data_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
		flag = system2(command=path_to_louvain, args=c("-i", data_path, "-o", result_path, "-r", resolution));
		if (flag != 0) {
		   	stop("'runCluster' call failed");
		}	
		cluster = read.table(result_path);
		object@cluster = factor(cluster[order(cluster[,1]),2]);	
		file.remove(data_path)
		file.remove(result_path)
	}else if(method=="kmeans"){
		if(is.null(centers)){
			stop("'centers' is null")
		}
		object@cluster = factor(kmeans(object@smat[,pca_dims], centers=centers)$cluster);
	}else if(method=="densityClust"){
		if(is.null(centers)){
			stop("'centers' is null")
		}
		irisDist <- dist(object@smat[,pca_dims]);
		irisClust <- densityClust(irisDist, gaussian=TRUE);
		cutoff_idx = order(irisClust$delta * irisClust$rho, decreasing=TRUE)[centers]
		delta_cutoff = irisClust$delta[cutoff_idx];
		rho_cutoff = irisClust$rho[cutoff_idx];		
		irisClust <- findClusters(irisClust, rho=rho_cutoff-0.0001, delta=delta_cutoff-0.0001);
		object@cluster = factor(irisClust$cluster);
	}
	return(object);
}





