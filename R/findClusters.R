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
	method=c("louvain", "jaccard_louvain", "kmeans", "densityClust"), 
	resolution=1, 
	data_path=NULL, 
	result_path=NULL,
	path_to_snaptools=NULL
){		
	if (is.null(data_path)) {
		data_path <- tempfile(pattern='community_data_', fileext='.dat')
	}
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='community_result_', fileext='.dat')
	}
	
	method = match.arg(method);
	if(method %in% c("jaccard_louvain", "louvain")){
		if (is.null(path_to_snaptools)) {
			path_to_snaptools <- system2('which', 'snaptools', stdout=TRUE)
		}
		path_to_snaptools <- normalizePath(path_to_snaptools);
		if (!file_test('-x', path_to_snaptools)) {
			stop(path_to_snaptools, " does not exist or is not executable; check your path_to_snaptools parameter")
		}		
	}
	
	# 1. check if object is a snap object
	if(class(object) != "snap"){
		stop("'object' is not a snap object");
	}

	# 2. check if smat exists
	if(nrow(object@smat) == 0){
		stop("smat does not exist in object, run PCA first");
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
		g = graph_from_edgelist(as.matrix(edges), directed=FALSE);
		#g = graph_from_adjacency_matrix(ajm, mode="undirected", weighted=NULL);		
		#g = graph_from_adjacency_matrix(similarity(g, method = "jaccard"), mode="undirected", weighted=TRUE)			
		#adj=get.adjacency(g,attr='weight') 
		gm = Matrix(similarity(g, method = "jaccard"), sparse=TRUE);		
		edges = data.frame(i= gm@i+1, j=findInterval(seq(gm@x)-1,gm@p[-1])+1, value=gm@x);
		edges = edges[edges[,1] < edges[,2],];
		write.table(edges, file = data_path, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
		flag = system2(command=path_to_snaptools, args=c("louvain", "--edge-file", data_path, "--output-file", result_path, "--resolution", resolution));
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
		flag = system2(command=path_to_snaptools, args=c("louvain", "--edge-file", data_path, "--output-file", result_path, "--resolution", resolution));
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





