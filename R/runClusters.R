#' Find Clusters Using Louvain Algorithm
#'
#' Using the constructed knn graph returned by function runKNN, we next applied community finding algorithm Louvain 
#' to identify the ‘communities’ in the resulting graph which represents groups of
#' cells sharing similar accessibility profiles.
#' 
#' @param obj A snap object.
#' @param tmp.folder Directory to store temporary files.
#' @param louvain.lib Louvain implementation method to use ["R-igraph", "python-louvain"].
#' "R-igraph" uses "cluster_louvain" implemented by igraph package in R. "python-louvain" uses louvain method implemented in "python-louvain" package. 
#' "R-igraph" is signficantly faster but does not support different resolutions. "python-louvain" supports
#' multiple resolutions, but requires "snaptools" and runs significantly slower.
#' @param resolution A numeric value that indicates the resolution for louvain clustering [1]. This is effective only when path.to.snaptools is set.
#' @param path.to.snaptools Path to snaptools excutable file [NULL]. Install snaptools by 'pip install snaptools'.
#' @param seed.use Random seed [10]. 
#' 
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(obj=demo.sp, tmp.folder=tempdir(), mat="bmat");
#' demo.sp = runNormJaccard(obj=demo.sp, tmp.folder=tempdir());
#' demo.sp = runDimReduct(obj=demo.sp, pc.num=10, input.mat="jmat");
#' demo.sp = runKNN(obj=demo.sp, pca.dims=1:5, k=15, snn=TRUE, save.knn=FALSE);
#' demo.sp = runCluster(obj=demo.sp, tmp.folder=tempdir(), louvain.lib="R-igraph");
#' 
#' @return Returns a snap obj with the cluster stored in obj@cluster
#'
#' @importFrom igraph graph_from_adjacency_matrix graph_from_edgelist cluster_louvain E
#' @importFrom utils file_test write.table read.table
#' @importFrom methods as
#' @export

runCluster <- function(obj, tmp.folder, louvain.lib, resolution, path.to.snaptools, seed.use) {
  UseMethod("runCluster", obj);
}

#' @export
runCluster.default <- function(
	obj, 
	tmp.folder, 
	louvain.lib=c("R-igraph", "python-louvain"),
	resolution=1.0,
	path.to.snaptools=NULL,
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
	
	if(isKgraphEmpty(obj@graph)){
		stop("knn graph is empty, run 'runKNN' first")		
	}
				
	if(is.na(as.numeric(resolution))){
		stop("resolution must be numeric class!")
	}
	
	if(louvain.lib == "R-igraph"){
		data.use = getGraph(obj@graph);
		g = graph_from_adjacency_matrix(data.use, weighted=TRUE, mode="undirected");
		cat("Epoch: finding clusters using R-igraph\n", file = stderr())		
		set.seed(seed.use);
		cl = cluster_louvain(g);
		obj@cluster = factor(cl$membership);		
	}else if(louvain.lib == "python-louvain"){
		if((x=nrow(obj@graph@mat)) > 0L){
			data_path <- tempfile(pattern='community_data_', tmpdir=tmp.folder, fileext='.dat');
			adj = obj@graph@mat;
			i = adj@i + 1;
			j = findInterval(seq(adj@x)-1,adj@p[-1])+1; 
			w = adj@x;
			edgeList = data.frame(v1=i, v2=j, w=w);
			writeEdgeListToFile(edgeList, data_path);						
		}else{
			data_path = obj@graph@file;
		}
		cat("Epoch: finding clusters using python-louvain\n", file = stderr())		
		result_path <- tempfile(pattern='community_result_', tmpdir=tmp.folder, fileext='.dat')
		flag = system2(command=path.to.snaptools, args=c("louvain", "--edge-file", data_path, "--output-file", result_path, "--resolution", resolution));		
		if (flag != 0) {
		   	stop("'runCluster' call failed");
		}				
		cluster = read.table(result_path);		
		file.remove(result_path);
		obj@cluster = factor(cluster[order(cluster[,1]),2]);		
		file.remove(result_path);
	}else{
		stop("unrecognized louvain.lib option")
	}
	gc();
	return(obj);
}


