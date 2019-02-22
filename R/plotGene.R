#' Plot Gene Accessibility Levels
#' @export

plotGene <- function(obj, ...) {
  UseMethod("plotGene", obj);
}

#' @export
plotGene.default <- function(
	obj, 
	gene.sel, 
	method=c("tsne", "umap"), 
	binary=TRUE,
	p=0.1, 
	plot.row=4, 
	plot.col=4, 
	background=TRUE,
	background_rho=0.2,
	...
	){	
	
		if(missing(obj) || missing(gene.sel)){
			stop("'obj' or 'gene.sel' is missing")
		}
		
		if(class(obj) != "snap"){
			stop("'obj' is not a snap object")
		}
		
		if(nrow(obj@gmat) == 0){
			stop("cell-by-gene matrix is empty")		
		}
		
		method = match.arg(method);
		if(method=="tsne"){
			if(nrow(obj@tsne) == 0){
				stop("tsne is missing")
			}
			x = obj@tsne;
		}

		if(method=="umap"){
			if(nrow(obj@umap) == 0){
				stop("umap is missing")
			}
			x = obj@umap;
		}
	
		if(any(!(gene.sel %in% colnames(obj@gmat)))){
			stop(paste(gene.sel[which(!(gene.sel %in% colnames(obj@gmat)))], "does not exist in cell-by-gene matrix"))
		};
		
		op <- par(mfrow = c(plot.row,plot.col), oma = c(3,3,1,1) + 0.2, mar = c(0,0,1,1) + 0.2);
		ncell = nrow(obj@gmat);			  
		ndim = length(gene.sel);
    	
		if(binary){
			for(i in seq(ndim)){
				id.sel = which(colnames(obj@gmat) == gene.sel[i]);
				y = obj@gmat[,id.sel];
				id.pos = order(y, decreasing=TRUE)[1:round(min(length(which(y > 0)), p*ncell))];	
				plot(x,
			   		col="grey", 
			   		main=gene.sel[i],
			   		yaxt='n', 
			   		xaxt="n",
			   		xlab="", 
			   		ylab="",
					...
			   		);
				points(x[id.pos,], col="red",...);
			}			
		}else{
			for(i in seq(ndim)){
				id.sel = which(colnames(obj@gmat) == gene.sel[i]);
				y = obj@gmat[,id.sel];
				if(background){
					plot(x, 
			   			 main=gene.sel[i],
						 col=alpha("grey", background_rho), 
	 			   		 yaxt='n', 
	 			   		 xaxt="n",
	 			   		 xlab="", 
	 			   		 ylab="",
	 					 ...
						 );					
					
					points(x, 
						 col=alpha("red", pmin(1, y/quantile(y[which(y > 0)], 0.99))), 
	 					 ...
						 );										
				}else{
					plot(x, 
			   			 main=gene.sel[i],
						 col=alpha("red", pmin(1, y/quantile(y[which(y > 0)], 0.99))), 
	 			   		 yaxt='n', 
	 			   		 xaxt="n",
	 			   		 xlab="", 
	 			   		 ylab="",
	 					 ...
						 );					
				}
				
			}						
		}
}

