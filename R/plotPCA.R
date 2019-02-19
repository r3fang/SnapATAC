#' PCA Plot
#'
#' @param obj a snap object
#' @export
plotPCA <- function(obj, ...) {
  UseMethod("plotPCA", obj);
}

#' @export
plotPCA.default <- function(object, method=c("elbow", "pairwise"), max.cell=5000){
	n <- 60
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
	col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("'object' is not a snap object")
	}

	# 2. check of jmat exists;
	if(nrow(object@smat) == 0){
		stop("@smat is empty")		
	}	
	
	method = match.arg(method);
	
	if(nrow(object@smat) > max.cell){
		idx = sort(sample(seq(nrow(object@smat)), max.cell))
	}else{
		idx = seq(nrow(object@smat));
	}
	
	if(method == "elbow"){
		plot(object@d,
	 	 	 xlab="PCs",
		 	 ylab="variance",
		 );
	}else{
		op <- par(mfrow = c(5,5), oma = c(3,3,1,1) + 0.2, mar = c(0,0,1,1) + 0.2)			  
		ndim = ncol(object@smat);
		if(length(object@cluster) == length(object@barcode)){
			for(i in seq(1, ndim, by=2)){
				plot(object@smat[idx,i], 
					 object@smat[idx,i+1], 
					 cex=0.1, 
					 col=col_vector[object@cluster][idx], 
 					 mtext(paste(paste("PC", i), i+1, sep=" vs "), side=3),
					 yaxt='n', 
					 xaxt="n",
					 xlab="", 
					 ylab=""
					 );
			}		
		}else{
			for(i in seq(1, ndim, by=2)){
				plot(object@smat[idx,i], 
					object@smat[idx,i+1], 
					cex=0.1, 
					col="grey",
					mtext(paste(paste("PC", i), i+1, sep=" vs "), side=3),
				 	yaxt='n', 
				 	xaxt="n",
				 	xlab="", 
				 	ylab=""
					);
			}				
		}		
	}
}

