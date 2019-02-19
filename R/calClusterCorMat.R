##' Calculate Cluster's Pearson Correlation Between Replicates and Clusters
##'
##' @param object A Snap object
##' @param cluster A numeric class that indicates the cluster label for each cell
##' @param rep A numeric class that indicates the replicate label for each cell
##'
##' @return Returns a numberic matrix contains pearwise correlation between any two cluster and replicates combination
##'
##' @export

#calClusterCorMat <- function(object, cluster, rep, ...) {
#  UseMethod("calClusterCorMat");
#}
#
#calClusterCorMat.default <- function(object, cluster=NULL, rep=NULL){
#	
#	# 1. check if object is a snap;
#	if(class(object) != "snap"){
#		stop("Not a snap object")
#	}
#		
#	if(missing(cluster)){
#		stop("cluster is missing")		
#	}
#
#	if(missing(rep)){
#		stop("rep is missing")		
#	}
#
#	if((nrow(object) != length(cluster)) | nrow(object) != length(rep)){
#		stop("cluster or rep incorrect length")				
#	}
#	cluster = factor(cluster)
#	rep = factor(rep)
#	
#	cluster = droplevels(cluster)
#	rep = droplevels(rep)
#	ncell = nrow(object)
#
#	# check if all cluster contains replicates
#	cluster_rep.df <- data.frame(cluster, rep, idx=1:ncell)
#	cluster_rep.num <- do.call(c, lapply(split(cluster_rep.df, cluster_rep.df$cluster), function(x){
#		x = droplevels(x)
#		length(levels(x[,2]))
#	}))
#	
#	# calculate average pearson corrlation between replicates for each cluster
#	avg.ls <- lapply(split(cluster_rep.df, cluster_rep.df$cluster), function(x){
#		x = droplevels(x)
#		levels.arr <- levels(x[,2])
#		if(length(levels.arr) < 2){
#			return(Matrix::colSums(object@cmat[x[,3],]))
#		}else{
#			cluster_rep_avg = matrix(0, ncol=length(levels.arr), nrow=ncol(object));
#			for(i in 1:length(levels.arr)){
#				cluster_rep_avg[,i] <- Matrix::colSums(object@cmat[x[which(x[,2] == levels.arr[i]),3],])
#			}
#			return(cluster_rep_avg)
#		}
#	})
#	
#	avg.mat <- do.call(cbind, avg.ls);
#	return(cor(avg.mat))
#}
#