#' Tsne Plot
#'
#' @param object A snap object.
#' @export

plotViz <- function(object, method, ...) {
  UseMethod("plotViz", object);
}

#' @export
plotViz.default <- function(object, method=c("tsne", "umap"), ...){
    tsnecols = c("#E31A1C","#FFD700","#771122","#777711","#1F78B4","#68228B","#AAAA44",
                 "#60CC52","#771155","#DDDD77","#774411","#AA7744","#AA4455","#117744",
                 "#000080","#44AA77","#AA4488","#DDAA77", colVector(100))
	
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("object is not a snap object")
	}
	
	method = match.arg(method);
	if(method=="tsne"){
		if(nrow(object@tsne) == 0L){
			stop("tsne does not exist, run runViz first!")
		}
		if(length(object@cluster) != 0L){
			plot(object@tsne, col=tsnecols[object@cluster], ..., main="tsne scatter plot", xlab="tsne-1", ylab="tsne-2")
			plot(object@tsne, col=tsnecols[object@cluster], ..., type="n", main="tsne cluster label", xlab="tsne-1", ylab="tsne-2")
			text(object@tsne, col=tsnecols[object@cluster], label=object@cluster, ...)
		}else{
			plot(object@tsne, col="grey", ...);
		}
	}else{
		if(nrow(object@umap) == 0L){
			stop("umap does not exist, run runViz first!")
		}
		if(length(object@cluster) != 0L){
			plot(object@umap, col=tsnecols[object@cluster], main="umap scatter plot", xlab="umap-1", ylab="umap-2", ...)
			plot(object@umap, col=tsnecols[object@cluster], type="n", main="umap cluster label", xlab="umap-1", ylab="umap-2", ...)
			text(object@umap, col=tsnecols[object@cluster], label=object@cluster, ...)
		}else{
			plot(object@umap, col="grey", ...);
		}		
	}
}

