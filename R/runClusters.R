#' Find Clusters Using Louvain Algorithm
#'
#' Using the constructed knn graph, we next applied community finding algorithm Louvain 
#' to identify the ‘communities’ in the resulting graph which represents groups of
#' cells sharing similar profiles, potentially originating from the same cell type.
#' 
#' @param obj A snap object
#' @param louvain.lib Louvain implementation method to use ["R-igraph", "python-louvain"].
#' "R-igraph" uses "cluster_louvain" to identify cluster. "python-louvain" uses "python-louvain" package to identify
#' clusters. "R-igraph" is signficantly faster but does not support different resolutions. "python-louvain" supports
#' multiple resolutions, but requires "snaptools" to be pre-installed.
#' @param resolution A numeric class that indicates the resolution for louvain clustering [1]. This is effective only when path.to.snaptools is set.
#' @param path.to.snaptools Path to snaptools excutable file [NULL]. Install snaptools by 'pip install snaptools'.
#' @param load.knn.from.file Whether to load knn grph from a file [FALSE].
#' @param edge.file Directory path to the file that contains knn graphs [NULL]. 
#' @param seed.use Random seed [10]. 
#' 
#' @return Returns a snap obj with the cluster stored in obj@cluster
#'
#' @importFrom igraph graph_from_adjacency_matrix graph_from_edgelist cluster_louvain E
#' @importFrom utils file_test write.table read.table
#' @importFrom methods as
#' @export

runCluster <- function(obj, louvain.lib, resolution, path.to.snaptools, load.knn.from.file, edge.file, seed.use) {
  UseMethod("runCluster", obj);
}

#' @export
runCluster.default <- function(
	obj, 
	louvain.lib=c("R-igraph", "python-louvain"),
	resolution=1.0,
	path.to.snaptools=NULL,
	load.knn.from.file=FALSE, 
	edge.file=NULL,
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

	louvain.lib = match.arg(louvain.lib);
	# louvain.lib is python-louvain only if 
	# 1) louvain.lib == "python-louvain";
	# 2) path.to.snaptools exists;
	# 3) path.to.snaptools excutable;
	if(louvain.lib == "python-louvain"){
		flag1 = (louvain.lib == "python-louvain");
		flag2 = tryCatch({
			file.exists(path.to.snaptools);		
		},
		error=function(cond){
			warning("path.to.snaptools does not exist, switch louvain.lib to R-igraph\n")
			return(FALSE)
		})
	
		flag3 = tryCatch({
			file_test('-x', path.to.snaptools);	
		},
		error=function(cond){
			warning("path.to.snaptools is not an excutable file, switch louvain.lib to R-igraph")
			return(FALSE)
		})
	
		if(flag1 && flag2 && flag3){
			louvain.lib = "python-louvain"
		}else{
			louvain.lib = "R-igraph";						
		}
	}
	
	if(nrow(obj@kmat) == 0 && load.knn.from.file==FALSE){
		stop("knn graph is empty, run 'runKNN' first")
	}
	
	if(nrow(obj@kmat) != 0 && load.knn.from.file==TRUE){
		warning("knn already exists in obj, will not read from edge.file")
		load.knn.from.file = FALSE;
	}

	if(load.knn.from.file==TRUE && is.null(edge.file)){
		stop("edge.file does not exist")
	}
			
	if(is.na(as.numeric(resolution))){
		stop("resolution must be numeric class!")
	}
	
	if(louvain.lib == "R-igraph"){
		if(load.knn.from.file){
			cat("Epoch: loading KNN graph from file\n", file = stderr())					
			if(!file.exists(edge.file)){
				stop("edge.file does not exist!")
			}
			edgeList = readEdgeListFromFile(edge.file);
			g = igraph::graph_from_edgelist(as.matrix(edgeList[,c(1,2)]), directed=FALSE);
			igraph::E(g)$weight = edgeList[,3]			
		}else{
			g = graph_from_adjacency_matrix(obj@kmat, weighted=TRUE, mode="undirected");
		}
		cat("Epoch: finding clusters using R-igraph\n", file = stderr())		
		set.seed(seed.use);
		cl = cluster_louvain(g);
		obj@cluster = factor(cl$membership);		
	}else if(louvain.lib == "python-louvain"){
		if(load.knn.from.file){
			if(!file.exists(edge.file)){
				stop("edge.file does not exist!")
			}
			data_path = edge.file;
		}else{
			data_path <- tempfile(pattern='community_data_', fileext='.dat');
			adj = obj@kmat;
			i = adj@i + 1;
			j = findInterval(seq(adj@x)-1,adj@p[-1])+1; 
			w = adj@x;
			edgeList = data.frame(v1=i, v2=j, w=w);
			writeEdgeListToFile(edgeList, data_path);			
		}
		cat("Epoch: finding clusters using python-louvain\n", file = stderr())		
		result_path <- tempfile(pattern='community_result_', fileext='.dat')
		flag = system2(command=path.to.snaptools, args=c("louvain", "--edge-file", data_path, "--output-file", result_path, "--resolution", resolution));		
		if (flag != 0) {
		   	stop("'runCluster' call failed");
		}				
		cluster = read.table(result_path);		
		file.remove(result_path);
		obj@cluster = factor(cluster[order(cluster[,1]),2]);		
	}else{
		stop("unrecognized louvain.lib option")
	}
	gc();
	return(obj);
}


