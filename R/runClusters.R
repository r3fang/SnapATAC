#' Find Clusters Using Louvain/Leiden Algorithm
#'
#' Using the constructed knn graph returned by function runKNN, we next applied community finding algorithm 
#' to identify the ‘communities’ in the resulting graph which represents groups of
#' cells sharing similar accessibility profiles.
#' 
#' @param obj A snap object.
#' @param tmp.folder Directory to store temporary files.
#' @param louvain.lib Louvain implementation method to use ["R-igraph", "leiden"].
#' "R-igraph" uses "cluster_louvain" implemented by igraph package in R. "Leiden" uses "Leiden" algorithm for finding clusters (recommanded). 
#' Leiden allows for multiple resolutions, but requires "leiden" to be pre-installed seperately. see how to install "leiden" (https://github.com/TomKellyGenetics/leiden).
#' @param resolution A numeric value that indicates the resolution for louvain clustering [1].
#' @param seed.use Random seed [10]. 
#' 
#' @examples
#' data(demo.sp);
#' demo.sp = runCluster(obj=demo.sp, tmp.folder=tempdir(), louvain.lib="R-igraph");
#' 
#' @return Returns a snap obj with the cluster stored in obj@cluster
#'
#' @importFrom igraph graph_from_adjacency_matrix graph_from_edgelist cluster_louvain E
#' @importFrom utils file_test write.table read.table
#' @importFrom methods as
#' @export
runCluster <- function(obj, tmp.folder, louvain.lib, resolution, seed.use) {
  UseMethod("runCluster", obj);
}

#' @export
runCluster.default <- function(
	obj, 
	tmp.folder, 
	louvain.lib=c("R-igraph", "leiden"),
	resolution=1.0,
	seed.use=10
){
	cat("Epoch: checking input parameters\n", file = stderr())
	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
	}
	
	if(missing(tmp.folder)){
		stop("tmp.folder is missing")
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist");			
		}
	}
	
	louvain.lib = match.arg(louvain.lib);
	
	if(louvain.lib == "leiden"){
		if (!requireNamespace("leiden", quietly = TRUE)) {
		      stop("Please install leiden - learn more at https://github.com/TomKellyGenetics/leiden")
		}
	}
	
	if(isKgraphEmpty(obj@graph)){
		stop("knn graph is empty, run 'runKNN' first")		
	}
				
	if(is.na(as.numeric(resolution))){
		stop("resolution must be numeric class!")
	}
	
	if(louvain.lib == "R-igraph"){
		data.use = getGraph(obj@graph);
		data.use = data.use + t(data.use);
		g = graph_from_adjacency_matrix(data.use, weighted=TRUE, mode="undirected");
		cat("Epoch: finding clusters using R-igraph\n", file = stderr())		
		set.seed(seed.use);
		cl = cluster_louvain(g);
		obj@cluster = factor(cl$membership);		
	}else if(louvain.lib == "leiden"){
		cat("Epoch: finding clusters using leiden\n", file = stderr())		
		data.use = getGraph(obj@graph);		
		set.seed(seed.use);
		obj@cluster <- factor(leiden(data.use, resolution=resolution));	
	}else{
		stop("unrecognized louvain.lib option")
	}
	gc();
	return(obj);
}

