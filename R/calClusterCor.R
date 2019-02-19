##' Calculate Cluster's Pearson Correlation Between Replicates
##'
##' One of the evidence that identified clusters are real is whether this cluster
##' is reproducible between biological replicates. This function is to calculate 
##' the pearson correlation between replicates for identified clusters
##' 
##' This function takes a snap object and a cluster vector and replicate vector as input
##' can calculate the pearon correlation between replicates for each cluster. The aeverage
##' correlation is returned when more than 2 replicates are given.
##' 
##' @param object A Snap object
##' @param cluster A numeric class that indicates the cluster label for each cell
##' @param rep A numeric class that indicates the replicate label for each cell
##'
##' @return Returns a numberic vector contains the correlation for each cluster
##'
##' @export

#calClusterCor <- function(object, cluster, rep, ...) {
#  UseMethod("calClusterCor");
#}
#
#calClusterCor.default <- function(object, cluster=NULL, rep=NULL){
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
#
#	ncell = nrow(object)
#	# check if all cluster contains replicates
#	cluster_rep.df <- data.frame(cluster, rep, idx=1:ncell)
#	cluster_rep.num <- do.call(c, lapply(split(cluster_rep.df, cluster_rep.df$cluster), function(x){
#		length(levels(x[,2]))
#	}))
#	
#	# if any cluster contains less than 2 cluster
#	if(any(cluster_rep.num < 2)){
#		message("Warning: found clusters does not contain replicates")
#	}
#	
#	# calculate average pearson corrlation between replicates for each cluster
#	cor.ls <- lapply(split(cluster_rep.df, cluster_rep.df$cluster), function(x){
#		levels.arr <- levels(x[,2])
#		if(length(levels.arr) < 2){
#			return(NA)
#		}else{
#			cluster_rep_avg = matrix(0, ncol=length(levels.arr), nrow=ncol(x.sel.sp));
#			for(i in 1:length(levels.arr)){
#				cluster_rep_avg[,i] <- Matrix::colSums(object@cmat[x[which(x[,2] == levels.arr[i]),3],])
#			}
#			cluster_rep_avg_cor = cor(cluster_rep_avg)
#			mean(cluster_rep_avg_cor[upper.tri(cluster_rep_avg_cor)])			
#		}
#	})
#	cor.arr <- do.call(c, cor.ls);
#	return(cor.arr)
#}
#