#' An S4 class snap to represent single-nucleus accessibility object.
#'
#' Class snap defines a snap object.
#'
#' @slot barcode A character vector contains cell barcodes as rows
#' @slot metaData A data.frame object contains meta data for barcodes
#' @slot feature A GRanges object contains genomic features (bins)
#' @slot peak A GRanges object contains genomic features (peaks)
#' @slot bmat A Matrix object contains cellxbin count matrix
#' @slot pmat A Matrix object contains cellxpeak count matrix
#' @slot gmat A matrix object contains cellxgene count matrix
#' @slot jmat A matrix object contains jaccard matrix
#' @slot nmat A matrix object contains normalized jaccard matrix
#' @slot smat A matrix object contains SVD projection matrix
#' @slot tsne A matrix object contains tsne coordinate
#' @slot umap A matrix object contains umap coordinate
#' @slot cluster A factor object contains cluster label
#' @name snap-class
#' @rdname snap-class
#' @exportClass snap

setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
setClass("snap",
	slots=list(
	barcode="character",
	metaData="data.frame",
	feature="GRanges",
	peak="GRanges",
	bmat = "Matrix",
	pmat = "Matrix",
	gmat = "Matrix",
	jmat = "MatrixOrmatrix",
	nmat = "MatrixOrmatrix",
	smat = "MatrixOrmatrix",
	graph = "MatrixOrmatrix",
	tsne = "MatrixOrmatrix",
	umap = "MatrixOrmatrix",
	cluster = "factor",
	d = "numeric"
	)
)

.valid.snap.feature <- function(object)
{
	if(length(object@feature) != ncol(object@bmat)){
		return("slot 'feature' have different length from 'bmat' column")		
	}
	NULL;
}

.valid.snap.peak <- function(object)
{
	if(length(object@peak) != ncol(object@pmat)){
		return("slot 'peak' have different length from 'pmat' column")		
	}
	NULL;
}

.valid.snap.barcode <- function(object)
{
	if(length(object@barcode) != nrow(object@metaData)){
		return("slot 'barcode' have different length from 'metaData'")		
	}
	NULL
}

.valid.snap <- function(object)
{
    #c(.valid.snap.barcode(object), .valid.snap.feature(object))
    c(.valid.snap.barcode(object), .valid.snap.peak(object), .valid.snap.feature(object));
}
setValidity("snap", .valid.snap)

#' Create a snap object
#' 
#' This function takes a matrix of single-nucleus count matrix 
#' (cell in rows and genomic features as columns) and a vector 
#' of barcodes and a GRanges of genomic features (genes or bins) 
#' to create a snap object
#'
#' @param cmat Matrix; A count matrix with barcodes as rows and genomic features as columns
#' @param barcode character; A vector object contains barcodes for the rows of mat
#' @param feature GRanges; A GRanges object contains genomic bins for the column of mat
#' @return A snap object
#'
#' @examples
#' barcode_num = 100
#' feature_num = 500
#' m <- matrix(sample(0:1,barcode_num * feature_num, replace=TRUE),barcode_num,feature_num)
#' newSnap(m)
#' 
#' @export
#' 
newSnap <- function (...) {
  UseMethod("newSnap");
}
newSnap.default <- function(){
	metaData=data.frame();
	barcode = as.character(c());
	feature = GRanges();
	peak = GRanges();		
	metaData = data.frame();
	bmat=Matrix(nrow=0, ncol=0, sparse=TRUE);		
	pmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	gmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	jmat=matrix(nrow=0, ncol=0);
	nmat=matrix(nrow=0, ncol=0);
	smat=matrix(nrow=0, ncol=0);	
	tsne=matrix(nrow=0, ncol=0);	
	umap=matrix(nrow=0, ncol=0);	
	graph=Matrix(nrow=0, ncol=0, sparse=TRUE);
	cluster=factor();
	res = new("snap", barcode=barcode, feature=feature, peak=peak, metaData=metaData, bmat=bmat, pmat=pmat, gmat=gmat, jmat=jmat, nmat=nmat, smat=smat, tsne=tsne, umap=umap, graph=graph, cluster=cluster, d=numeric(0));		
}

#newSnap.default <- function(cmat=NULL, barcode=NULL, feature=NULL, peak=NULL, metaData=NULL, bmat=NULL, pmat=NULL, gmat=NULL, jmat=NULL, nmat=NULL, smat=NULL, tsne=NULL, umap=NULL, graph=NULL, cluster=NULL, idx=NULL, cov=NULL){
#	# 
#
#	# empty initlization
#	if(missing(cmat) && missing(barcode) && missing(feature)){
#		metaData=data.frame();
#		barcode = as.character(c());
#		feature = GRanges()
#		metaData = data.frame();
#		cmat=Matrix(nrow=0, ncol=0, sparse=TRUE);		
#		bmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
#		jmat=matrix(nrow=0, ncol=0);
#		nmat=matrix(nrow=0, ncol=0);
#		gmat=matrix(nrow=0, ncol=0);
#		pmat=matrix(nrow=0, ncol=0);	
#		tsne=matrix(nrow=0, ncol=0);	
#		umap=matrix(nrow=0, ncol=0);	
#		graph=Matrix(nrow=0, ncol=0, sparse=TRUE);
#		cluster=factor();
#		cov=factor();
#		idx=numeric();
#		res = new("snap", cmat=cmat, barcode=barcode, feature=feature, metaData=metaData, bmat=bmat, jmat=jmat, nmat=nmat, pmat=pmat, tsne=tsne, umap=umap, graph=graph, cluster=cluster, gmat=gmat, d=numeric(0), idx=idx, cov=cov);
#	}else{
#		# cmat, barcode and feature are required to create a snap object
#		if(missing(cmat)){stop("Error @newSnap: cmat is missing!")};
#		if(missing(barcode)){stop("Error @newSnap: barcode is missing!")};
#		if(missing(feature)){stop("Error @newSnap: feature is missing!")};
#
#		# check input format
#		if(!(class(cmat) %in% c("dgCMatrix", 'dgeMatrix', "matrix"))){stop("Error @snap: cmat is not a matrix!")}
#		if(class(cmat) != "sparseMatrix"){cmat = as(cmat, "sparseMatrix")};				
#	
#		# check validity of barcode
#		if(class(barcode) != "character"){stop("Error @newSnap: barcode is not character object!")};
#		if(length(barcode) != nrow(cmat)){stop("Error @newSnap: barcode length does not match with count matrix row number!")};
#		if(length(unique(barcode)) != length(barcode)){stop("Error @newSnap: duplicate barcode identified!")};
#
#		if(class(feature) != "GRanges"){stop("Error @snap: feature is not GRanges object!")}
#		if(length(feature) != ncol(cmat)){stop("Error @newSnap: feature length does not match with count matrix column number!")}
#
#		if(is.null(metaData)){metaData=data.frame()};
#		if(is.null(peak)){peak=GRanges()};	
#		if(is.null(bmat)){bmat=Matrix(nrow=0, ncol=0, sparse=TRUE)};
#		if(is.null(jmat)){jmat=matrix(nrow=0, ncol=0)};
#		if(is.null(nmat)){nmat=matrix(nrow=0, ncol=0)};
#		if(is.null(gmat)){gmat=matrix(nrow=0, ncol=0)};
#		if(is.null(pmat)){pmat=matrix(nrow=0, ncol=0)};	
#		if(is.null(tsne)){tsne=matrix(nrow=0, ncol=0)};	
#		if(is.null(umap)){umap=matrix(nrow=0, ncol=0)};	
#		if(is.null(graph)){graph=Matrix(nrow=0, ncol=0, sparse=TRUE)};
#		if(is.null(cluster)){cluster=factor()};
#		if(is.null(cov)){cov=factor()};
#		if(is.null(idx)){
#			if(nrow(cmat) == 0){
#				idx = numeric();
#			}else{
#				idx = 1:nrow(cmat);
#			}
#		};
#		res = new("snap", cmat=cmat, barcode=barcode, feature=feature, peak=peak, metaData=metaData, bmat=bmat, jmat=jmat, nmat=nmat, pmat=pmat, tsne=tsne, umap=umap, graph=graph, cluster=cluster, gmat=gmat, d=numeric(0), idx=idx, cov=cov);
#	}
#}

#' getBarcode for snap object.
#'
#' This function takes a snap object and returns the barcode information
#' @name getBarcode
#' @param obj A snap object
#' @examples
#' fname <- system.file("extdata", "demo_snap", package = "SNAPATAC")
#' load(fname)
#' barcodes <- getBarcode(schep.sp)
#' @rdname getBarcode-methods
#' @exportMethod getBarcode
setGeneric("getBarcode", function(obj) standardGeneric("getBarcode"))

#' @rdname getBarcode-methods
#' @aliases getBarcode,snap-method
setMethod("getBarcode", "snap", function(obj) obj@barcode)

#' getMetaData for snap object.
#'
#' This function takes a snap object and returns the metaData information
#' @name getMetaData
#' @param obj A snap object
#' @examples
#' fname <- system.file("extdata", "demo_snap", package = "SNAPATAC")
#' load(fname)
#' barcodes <- getMetaData(schep.sp)
#' @rdname getMetaData-methods
#' @exportMethod getMetaData
setGeneric("getMetaData", function(obj) standardGeneric("getMetaData"))

#' @rdname getBarcode-methods
#' @aliases getBarcode,snap-method
setMethod("getMetaData", "snap", function(obj) obj@metaData)

#' getFeature for snap object.
#'
#' This function takes a snap object and returns the features
#' @name getFeature
#' @param obj snap; a snap object
#' @examples
#' barcode_num <- 100
#' feature_num <- 500
#' m <- matrix(sample(0:1,barcode_num * feature_num, replace=TRUE),barcode_num,feature_num)
#' x.sp <- newSnap(m)
#' features <- getFeature(x.sp)
#' @rdname getFeature-methods
#' @exportMethod getFeature
setGeneric("getFeature", function(obj) standardGeneric("getFeature"))

#' @rdname getFeature-methods
#' @aliases getFeature,snap-method
setMethod("getFeature", "snap", function(obj) obj@feature)

#' getBinMatrix for snap object.
#'
#' This function takes a snap object and returns the cellxbin matrix
#' @name getBinMatrix
#' @param obj snap; a snap object
#' @rdname getBinMatrix-methods
#' @exportMethod getBinMatrix
setGeneric("getBinMatrix", function(obj) standardGeneric("getBinMatrix"))

#' @rdname getBinMatrix-methods
#' @aliases getBinMatrix,snap-method
setMethod("getBinMatrix", "snap", function(obj) obj@bmat)

#' getPeakMatrix for snap object.
#'
#' This function takes a snap object and returns the cellxpeak matrix
#' @name getPeakMatrix
#' @param obj snap; a snap object
#' @rdname getPeakMatrix-methods
#' @exportMethod getPeakMatrix
setGeneric("getPeakMatrix", function(obj) standardGeneric("getPeakMatrix"))

#' @rdname getPeakMatrix-methods
#' @aliases getPeakMatrix,snap-method
setMethod("getPeakMatrix", "snap", function(obj) obj@pmat)

#' getJaccardMatrix for snap object.
#'
#' This function takes a snap object and returns the jaccard matrix
#' @name getJaccardMatrix
#' @param obj snap; a snap object
#' @rdname getJaccardMatrix-methods
#' @exportMethod getJaccardMatrix
setGeneric("getJaccardMatrix", function(obj) standardGeneric("getJaccardMatrix"))

#' @rdname getJaccardMatrix-methods
#' @aliases getJaccardMatrix,snap-method
setMethod("getJaccardMatrix", "snap", function(obj) obj@jmat)

#' getPCA for snap object.
#'
#' This function takes a snap object and returns the pca cell.embedding
#' @name getPCA
#' @param obj snap; a snap object
#' @rdname getPCA-methods
#' @exportMethod getPCA
setGeneric("getPCA", function(obj) standardGeneric("getPCA"))

#' @rdname getPCA-methods
#' @aliases getPCA,snap-method
setMethod("getPCA", "snap", function(obj) obj@smat)

#' getTSNE for snap object.
#'
#' This function takes a snap object and returns its tsne coordinate
#' @name getTSNE
#' @param obj snap; a snap object
#' @rdname getTSNE-methods
#' @exportMethod getTSNE
setGeneric("getTSNE", function(obj) standardGeneric("getTSNE"))

#' @rdname getTSNE-methods
#' @aliases getTSNE,snap-method
setMethod("getTSNE", "snap", function(obj) obj@tsne)

#' getUMAP for snap object.
#'
#' This function takes a snap object and returns its tsne coordinate
#' @name getUMAP
#' @param obj snap; a snap object
#' @rdname getUMAP-methods
#' @exportMethod getUMAP
setGeneric("getUMAP", function(obj) standardGeneric("getUMAP"))

#' @rdname getUMAP-methods
#' @aliases getUMAP,snap-method
setMethod("getUMAP", "snap", function(obj) obj@umap)

#' nrow for snap object.
#'
#' This function takes a snap object and returns number of barcodes
#' @name nrow
#' @param x snap; a snap object
#' @rdname nrow-methods
#' @aliases nrow,snap-method
#' @exportMethod nrow
setMethod("nrow", "snap", function(x) length(x@barcode));

#' ncol for snap object.
#'
#' This function takes a snap object and returns number of features
#' @name ncol
#' @param x snap; a snap object
#' @rdname ncol-methods
#' @aliases ncol,snap-method
#' @exportMethod ncol
setMethod("ncol", "snap", function(x) length(x@feature));

#' rowSums for snap objects
#'
#' This function takes a snap object and returns the row sums of its count matrix.
#' @name rowSums
#' @param x snap; a snap object
#' @rdname rowSums-methods
#' @aliases rowSums,snap-method
#' @exportMethod rowSums
setMethod("rowSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::rowSums(x@bmat, na.rm);		
	}else if(mat == "pmat"){
		res = Matrix::rowSums(x@pmat, na.rm);				
	}else if(mat == "gmat"){
		res = Matrix::rowSums(x@gmat, na.rm);
	}
	return(res);
});

#' colSums for snap objects
#'
#' This function takes a snap object and returns the column sums of its count matrix.
#' @name colSums
#' @param x snap; a snap object
#' @rdname colSums-methods
#' @aliases colSums,snap-method
#' @exportMethod rowSums
setMethod("colSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::colSums(x@bmat, na.rm);		
	}else if(mat == "pmat"){
		res = Matrix::colSums(x@pmat, na.rm);				
	}else if(mat == "gmat"){
		res = Matrix::colSums(x@gmat, na.rm);
	}
	return(res);
});

#' rowMeans for snap objects
#'
#' This function takes a snap object and returns the row means of its count matrix.
#' @name rowMeans
#' @param x snap; a snap object
#' @rdname rowMeans-methods
#' @aliases rowMeans,snap-method
#' @exportMethod rowMeans
setMethod("rowMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::rowMeans(x@bmat, na.rm);		
	}else if(mat == "pmat"){
		res = Matrix::rowMeans(x@pmat, na.rm);				
	}else if(mat == "gmat"){
		res = Matrix::rowMeans(x@gmat, na.rm);
	}
	return(res);
});


#' colMeans for snap objects
#'
#' This function takes a snap object and returns the column means of its count matrix.
#' @name colMeans
#' @param x snap; a snap object
#' @rdname colMeans-methods
#' @aliases colMeans,snap-method
#' @exportMethod colMeans
setMethod("colMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::colMeans(x@bmat, na.rm);		
	}else if(mat == "pmat"){
		res = Matrix::colMeans(x@pmat, na.rm);				
	}else if(mat == "gmat"){
		res = Matrix::colMeans(x@gmat, na.rm);
	}
	return(res);
});

#' rowCovs for snap object.
#'
#' This function takes a snap object and returns the number of non-zero elements for each row of the count matrix
#' @name rowCovs
#' @param obj snap; a snap object
#' @examples
#' barcode_num <- 100
#' feature_num <- 500
#' m <- matrix(sample(0:1,barcode_num * feature_num, replace=TRUE),barcode_num,feature_num)
#' x.sp <- newSnap(m)
#' barcode_cov <- rowCovs(x.sp)
#' @rdname rowCovs-methods
#' @exportMethod rowCovs
setGeneric("rowCovs", function(x, ...) standardGeneric("rowCovs"))

#' @rdname rowCovs-methods
#' @aliases rowCovs,snap-method
setMethod("rowCovs", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::rowSums(obj@bmat != 0, na.rm);
	}else if(mat == "pmat"){
		res = Matrix::rowSums(obj@pmat != 0, na.rm);
	}else if(mat == "gmat"){
		res = Matrix::rowSums(obj@gmat != 0, na.rm);
	}
	return(res);
});

#' colCovs for snap object.
#'
#' This function takes a snap object and returns the number of non-zero elements for each column of the count matrix
#' @name colCovs
#' @param obj snap; a snap object
#' @examples
#' barcode_num <- 100
#' feature_num <- 500
#' m <- matrix(sample(0:1,barcode_num * feature_num, replace=TRUE),barcode_num,feature_num)
#' x.sp <- newSnap(m)
#' feature_cov <- colCovs(x.sp)
#' @rdname colCovs-methods
#' @exportMethod colCovs
setGeneric("colCovs", function(x, ...) standardGeneric("colCovs"))

#' @rdname rowCovs-methods
#' @aliases rowCovs,snap-method
setMethod("colCovs", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=FALSE){
	mat = match.arg(mat);
	if(mat == "bmat"){
		res = Matrix::colSums(obj@bmat != 0, na.rm);
	}else if(mat == "pmat"){
		res = Matrix::colSums(obj@pmat != 0, na.rm);
	}else if(mat == "gmat"){
		res = Matrix::colSums(obj@gmat != 0, na.rm);
	}
	return(res);
});

#' isSnap for snap object.
#'
#' Functions to check if an object is a snap object
#' @name is.snap
#' @param x any R object.
#' @rdname is.snap-methods
#' @exportMethod is.snap
setGeneric("is.snap", function(x) standardGeneric("is.snap"))

#' @rdname is.snap-methods
#' @aliases is.snap,snap-method
setMethod("is.snap", "snap", function(x) return(class(x) == "snap"));

#' findOverlaps for snap objects
#'
#' This function takes a snap object and returns the column sums of its count matrix.
#' @name colSums
#' @param x snap; a snap object
#' @rdname colSums-methods
#' @aliases colSums,snap-method
#' @exportMethod rowSums
setGeneric("FindOverlaps", function(x, y, ...) standardGeneric("FindOverlaps"))
setMethod("FindOverlaps", "snap", function(x, y, mat=c("bmat", "pmat")){
	mat = match.arg(mat);
	if(mat == "bmat"){
		ov = findOverlaps(x@feature, y);
	}else if(mat == "pmat"){
		ov = findOverlaps(x@peak, y);
	}
	return(ov);
});

#' rmSlot for snap object.
#'
#' Function to remove a certain slot from a snap object
#' @name rmSlot
#' @param x any snap object.
#' @rdname rmSlot-methods
#' @exportMethod rmSlot
setGeneric("rmSlot", function(x, ...) standardGeneric("rmSlot"))

#' @rmSlot rmSlot-methods
#' @aliases rmSlot,snap-method
setMethod("rmSlot", "snap", function(x, mat=c("bmat", "pmat", "gmat", "jmat", "nmat", "smat", "tsne", "umap", "cluster")){
	mat = match.arg(mat);
	if(mat == "bmat"){
		x@bmat = Matrix(0, 0, 0);
	}else if(mat == "gmat"){
		x@gmat = Matrix(0, 0, 0);
	}else if(mat == "pmat"){
		x@pmat = Matrix(0, 0, 0);
	}else if(mat == "jmat"){
		x@jmat = matrix(0, 0, 0);
	}else if(mat == "nmat"){
		x@nmat = matrix(0, 0, 0);		
	}else if(mat == "smat"){
		x@smat = matrix(0, 0, 0);		
	}else if(mat == "tsne"){
		x@tsne = matrix(0, 0, 0);		
	}else if(mat == "umap"){
		x@umap = matrix(0, 0, 0);
	}else if(mat == "cluster"){
		x@cluster = c();
	}
	return(x);
});

setMethod("show", signature = "snap",
	definition = function(object) {
		cat("number of barcodes: ", ifelse(is.null(length(object@barcode)), 0, length(object@barcode)), "\n", sep="");
		cat("number of bins: ", ncol(object@bmat), "\n", sep="");
		cat("number of peaks: ", ncol(object@pmat), "\n", sep="");
		cat("number of genes: ", ncol(object@gmat), "\n", sep="");
		cat("==========================\n");
		cat("meta data            (metData) : ", nrow(object@metaData) != 0, "\n");
		cat("cellxbin matrix      (bmat)    : ", nrow(object@bmat) != 0, "\n");
		cat("cellxpeak matrix     (pmat)    : ", nrow(object@pmat) != 0, "\n");
		cat("cellxgene matrix     (gmat)    : ", nrow(object@gmat) != 0, "\n");
		cat("jaccard matrix       (jmat)    : ", nrow(object@jmat) != 0, "\n");
		cat("normalization        (nmat)    : ", nrow(object@nmat) != 0, "\n");
		cat("PCA:                 (smat)    : ", nrow(object@smat) != 0, "\n");
		cat("cluster:             (cluster) : ", length(object@cluster) != 0, "\n");
		cat("t-sne:               (tsne)    : ", nrow(object@tsne) != 0, "\n");
		cat("umap:                (umap)    : ", nrow(object@umap) != 0, "\n");
	}
)

#' subsetting for snap objects
#'
#' This function takes a snap object and returns the subset of snap object.
#' @param x snap; a snap object
#' @param i integer; selected barcode index
#' @param j integer; selected feature index
#' @export
setMethod("[", "snap",
	function(x,i,j,mat=c("bmat", "pmat"), drop="missing"){
		.barcode = x@barcode;
		.feature = x@feature;
		.peak = x@peak;
		.bmat = x@bmat;
		.pmat = x@pmat;
		.gmat = x@gmat;
		.jmat = x@jmat;
		.nmat = x@nmat;
		.smat = x@smat;
		.cluster = x@cluster;
		.tsne = x@tsne;
		.umap = x@umap;
		.metaData = x@metaData;
		# a single row or column
       if(!missing(i)){
		   if(nrow(.bmat) > 0){.bmat <- .bmat[i,,drop=FALSE]}
		   if(nrow(.pmat) > 0){.pmat <- .pmat[i,,drop=FALSE]}
		   if(nrow(.gmat) > 0){.gmat <- .gmat[i,,drop=FALSE]}	   
		   if(nrow(.jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
		   if(nrow(.nmat) > 0){.nmat <- .nmat[i,,drop=FALSE]}
		   if(nrow(.smat) > 0){.smat <- .smat[i,,drop=FALSE]}
		   if(nrow(.tsne) > 0){.tsne <- .tsne[i,,drop=FALSE]}
		   if(nrow(.umap) > 0){.umap <- .umap[i,,drop=FALSE]}
		   if(nrow(.metaData) > 0){.metaData <- .metaData[i,,drop=FALSE]}
		   if(length(.cluster) > 0){.cluster <- .cluster[i,drop=FALSE]}
		   .barcode <- .barcode[i];
	   }
	   if(!missing(j)){
   			mat = match.arg(mat);
	   		if(mat == "bmat"){
	 		   if(nrow(.bmat) > 0){.bmat <- .bmat[,j,drop=FALSE]}
			   if(length(.feature) > 0){.feature <- .feature[j];}	   
	   		}else if(mat == "pmat"){
 	 		   if(nrow(.pmat) > 0){.pmat <- .pmat[,j,drop=FALSE]}
			   if(length(.peak) > 0){.peak <- .peak[j];}	   
	   		}
	   }
	   x@bmat = .bmat;
	   x@pmat = .pmat;
	   x@gmat = .gmat;
	   x@barcode = .barcode;
	   x@peak = .peak;
	   x@feature = .feature;
	   x@metaData = .metaData;
	   x@umap = .umap;
	   x@feature = .feature;
	   x@jmat = .jmat;
	   x@nmat = .nmat;
	   x@smat = .smat;
	   x@cluster = .cluster;
	   x@tsne = .tsne;
	   return(x);
})
