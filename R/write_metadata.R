#' Write down the meta data as text file
#' @export

write.metadata <- function(obj, file, ...){
    UseMethod("write.metadata", obj);	
}


#' @export
write.metadata.default <- function(obj, file, subset.names){
	subset.names.ref = c('barcode', 'cluster', 'tsne', 'umap', 'qc');
	
	if(class(obj) != "snap"){
		stop("'obj' is not a snap object")
	};
	
	if(file.exists(file)){
		stop("file already exists, remove it first")		
	}
	
	if(missing(subset.names)){
		stop("subset.names is missing")
	}

	if(!all(subset.names %in% subset.names.ref)){
		stop("'subset.names' must be subset of c('barcode', 'cluster', 'tsne', 'umap', 'qc')");		
	}
		
	metaData.ls = lapply(as.list(subset.names), function(x){
		if(x == "barcode"){
			y = data.frame(slot(obj, x));			
			colnames(y) = "barcode"
		}else if(x == "tsne"){
			y = data.frame(slot(obj, x));			
			colnames(y) = c("tsne1", "tsne2");
		}else if(x == "umap"){
			y = data.frame(slot(obj, x));			
			colnames(y) = c("umap1", "umap2");
		}else if(x == "cluster"){
			y = data.frame(slot(obj, x));			
			colnames(y) = "cluster"
		}else{
			y = data.frame(slot(obj, x));			
		}
		y
	})
	if(!all(sapply(lapply(metaData.ls, nrow), FUN = identical, nrow(metaData.ls[[1]])))){
		stop("slot in subset.names have different length")
	}
	
	metaData.df = do.call(cbind, metaData.ls);
	

    write.table(metaData.df, file = file, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = TRUE, qmethod = c("escape", "double"),
                fileEncoding = "")
}