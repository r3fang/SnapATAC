#' Convert a snap object to a seurat object
#'
#' @param obj A snap object.
#' @param eigs.dims A vector of the dimensions to use.
#' @param norm A logical variable indicates whether to normalize the cell by gene matrix.
#' @param scale A logical variable indicagtes whether to scale the cell by gene matrix.
#' @return Return a seurat object that contains single cell ATAC-seq data
#'
#' @export
snapToSeurat <- function(
	obj, 
	eigs.dims=1:20,
	norm=TRUE,
	scale=TRUE
){
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
		if((x=length(obj@barcode))==0L){
			stop("obj@barcode is empty");			
		}
		if((x=length(obj@file))==0L){
			stop("obj@file is empty");			
		}
	}
	# check if Seurat is installed
	if (requireNamespace("Seurat", quietly = TRUE)) {
		require(Seurat)
	} else {
		stop("Please install Seurat V3 - learn more at https://github.com/satijalab/seurat");
	}
	 
	 if((x=nrow(obj@gmat)) == 0L){
		 stop("gmat in obj is empty!");
	 }
	 gmat.use = t(obj@gmat);
 	# check if mat is binary;
	 
	if((x=nrow(obj@bmat)) > 0L){
		input.mat = "bmat";
	}else if((x=nrow(obj@pmat)) > 0L){
		input.mat = "pmat";
	}else{
		stop("both pmat and bmat is empty")
	}

	 if(input.mat == "bmat"){
		 data.use = obj@bmat;
		 peak.use = as.data.frame(obj@feature);
	 }else{
		 data.use = obj@pmat;
		 peak.use = as.data.frame(obj@peak);
	 }
	 
	 if((x=nrow(data.use)) == 0L){
		 stop("input matrix is empty!")
	 }
	 
	 metaData.use = obj@metaData;
	 if((x=nrow(metaData.use)) == 0L){
		 stop("metaData is empty!")	 	
	 }
	 
 	ncell = nrow(obj);
 	nvar = ncol(obj@smat@dmat);
	
 	if(missing(eigs.dims)){
 		stop("eigs.dims is missing")
 	}else{
 		if(is.null(eigs.dims)){
 			eigs.dims=1:nvar;	
 		}else{
 			if(any(eigs.dims > nvar) ){
 				stop("'eigs.dims' exceeds PCA dimentions number");
 			}		
 		}
 	}
	 
	 pca.use = obj@smat@dmat;
	 if((x=nrow(pca.use)) == 0L){
		 stop("dimentionality reduction is empty, runLDM first")
	 }else{
	 	pca.use = pca.use[,eigs.dims]
	 }
	
	 data.use = t(data.use);
	 rownames(x = data.use) = peak.use$name;
	 colnames(x = data.use) = paste0(obj@barcode, 1:ncell);
	 colnames(x = gmat.use) = paste0(obj@barcode, 1:ncell);
	 rownames(x = pca.use)  = paste0(obj@barcode, 1:ncell);
	 rownames(metaData.use) = paste0(obj@barcode, 1:ncell);
	 
	 pbmc.atac <- CreateSeuratObject(counts = data.use[1:10,], assay = "ATAC");
	 pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gmat.use);
	 pbmc.atac <- AddMetaData(pbmc.atac, metadata = metaData.use);
	 pbmc.atac$tech <- "atac"
	 DefaultAssay(pbmc.atac) <- "ATAC";
	 
	 colnames(x = pca.use) <- paste0("DC_", eigs.dims);
	 pbmc.atac[["SnapATAC"]] <- new(Class = "DimReduc", cell.embeddings = pca.use,
	       feature.loadings = matrix(0,0,0), feature.loadings.projected = matrix(0,0,0),
	       assay.used ="ATAC", stdev = rep(1,length(eigs.dims)), 
		   key ="DC_", jackstraw = new(Class = "JackStrawData"), misc = list()) 
	
	 DefaultAssay(pbmc.atac) <- "ACTIVITY";
	 if(norm){
		 pbmc.atac <- NormalizeData(pbmc.atac);	 	
	 }
	 if(scale){
		 pbmc.atac <- ScaleData(pbmc.atac);
	 }	
	return(pbmc.atac);
}
