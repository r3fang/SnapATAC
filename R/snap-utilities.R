#' @importFrom methods new
newSnap <- function () {
	metaData=data.frame();
	des = character()
	file = as.character(c());
	barcode = as.character(c());
	feature = GRanges();
	peak = GRanges();		
	metaData = data.frame();
	bmat=Matrix(nrow=0, ncol=0, sparse=TRUE);		
	pmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	gmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	jmat=newJaccard();
	smat=newDimReduct();	
	graph=newKgraph();
	tsne=matrix(nrow=0, ncol=0);	
	umap=matrix(nrow=0, ncol=0);	
	cluster=factor();
	res = new("snap", 
			  des=des,
			  file=file,
			  barcode=barcode, 
			  feature=feature, 
			  peak=peak, 
			  metaData=metaData, 
			  bmat=bmat, 
			  pmat=pmat, 
			  gmat=gmat, 
			  jmat=jmat, 
			  smat=smat, 
			  graph=graph, 
			  tsne=tsne, 
			  umap=umap, 
			  cluster=cluster
			  );	
}

#' summarySnap for snap object.
#'
#' This function takes a snap object and returns the summary statistics
#' @name summarySnap
#' @param obj snap; a snap object
#' @rdname summarySnap-methods
#' @importFrom stats median
#' @exportMethod summarySnap
setGeneric("summarySnap", function(obj) standardGeneric("summarySnap"))

#' @rdname summarySnap-methods
#' @aliases summarySnap,snap-method
setMethod("summarySnap", "snap", function(obj){
	if(nrow(obj@metaData) == 0){stop("metaData is empty")}
	barcode = obj@metaData;
	message("Total  number of barcodes: ", length(obj@barcode));	
	message("Median number of sequencing fragments: ", median(barcode$TN));
	message("Median number of uniquely mapped fragments: ", median(barcode$UQ));
	message("Median number of mappability ratio: ", round(median((barcode$UM+1)/(barcode$TN+1)),2));
	message("Median number of properly paired ratio: ", round(median((barcode$PP+1)/(barcode$UM+1)),2));
	message("Median number of duplicate ratio: ", round(median(1 - (barcode$UQ+1)/(barcode$PP+1)),2));
	message("Median number of chrM ratio: ", round(median((barcode$CM+1) / (barcode$UQ+1)),2));
	message("Median number of unique molecules (UMI): ", median(barcode$UQ));
});


#' nrow for snap object.
#'
#' This function takes a snap object and returns number of cells
#' @name nrow
#' @param x snap; a snap object
#' @rdname nrow-methods
#' @aliases nrow,snap-method
#' @exportMethod nrow
setMethod("nrow", "snap", function(x) length(x@barcode));

#' colSums for snap objects
#'
#' This function takes a snap object and returns the column sums of its count matrix.
#' @name colSums
#' @param x A snap object
#' @param mat A charater object indicates what matrix slot to use c("bmat", "pmat", "gmat")
#' @param na.rm A logical variable indicates wether to remove NA in the matrix
#' @rdname colSums-methods
#' @aliases colSums,snap-method
#' @exportMethod colSums
setMethod("colSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
	mat = match.arg(mat);
	mat.use = methods::slot(x, mat);
	if((x=nrow(mat.use))==0L){
		stop("mat is empty")
	}
	res = Matrix::colSums(mat.use, na.rm);
	return(res);
});

#' rowSums for snap objects
#'
#' This function takes a snap object and returns the row sums of its count matrix.
#' @name rowSums
#' @param x A snap object
#' @param mat A charater object indicates what matrix slot to use
#' @param na.rm A logical variable indicates wether to remove NA in the matrix
#' @rdname rowSums-methods
#' @aliases rowSums,snap-method
#' @exportMethod rowSums
setMethod("rowSums", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
	mat = match.arg(mat);
	mat.use = methods::slot(x, mat);
	if((x=nrow(mat.use))==0L){
		stop("mat is empty")
	}
	res = Matrix::rowSums(mat.use, na.rm);
	return(res);
});


#' rowMeans for snap objects
#'
#' This function takes a snap object and returns the row means of its count matrix.
#' 
#' @name rowMeans
#' @param x A snap object
#' @param mat A charater object indicates what matrix slot to use for calculation
#' @param na.rm A logical variable indicates wether to remove NA in the matrix
#' @rdname rowMeans-methods
#' @importFrom methods slot
#' @aliases rowMeans,snap-method
#' @exportMethod rowMeans
setMethod("rowMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
	mat = match.arg(mat);
	mat.use = methods::slot(x, mat);
	if((x=nrow(mat.use))==0L){
		stop("mat is empty")
	}
	res = Matrix::rowMeans(mat.use, na.rm);
	return(res);
});


#' colMeans for snap objects
#'
#' This function takes a snap object and returns the column means of its count matrix.
#' 
#' @name colMeans
#' @param x A snap object
#' @param mat A charater object indicates what matrix slot to use c("bmat", "pmat", "gmat").
#' @param na.rm A logical variable indicates wether to remove NA in the matrix.
#' @rdname colMeans-methods
#' @importFrom methods slot
#' @aliases colMeans,snap-method
#' @exportMethod colMeans
setMethod("colMeans", "snap", function(x, mat=c("bmat", "pmat", "gmat"), na.rm=TRUE){
	mat = match.arg(mat);
	mat.use = methods::slot(x, mat);
	if((x=nrow(mat.use))==0L){
		stop("mat is empty")
	}
	res = Matrix::colMeans(mat.use, na.rm);
	return(res);
});

#' Check snap object
#'
#' This function takes any object as input and check if it is a snap object
#' @param obj A snap object
#' @rdname is.snap-methods
#' @exportMethod is.snap
setGeneric("is.snap", function(obj) standardGeneric("is.snap"))

#' @rdname is.snap-methods
#' @aliases is.snap,snap-method
setMethod("is.snap", "snap", function(obj) return(class(obj) == "snap"));

#' subsetting for snap objects
#'
#' This function takes a snap object and returns the subset of snap object.
#' @param x snap; a snap object
#' @param i integer; selected barcode index
#' @param j integer; selected feature index
#' @param mat character; indicates the slot to subsetting
#' @param drop character; 
#' @export
setMethod("[", "snap",
	function(x,i,j,mat=c("bmat", "pmat"), drop="missing"){
		.barcode = x@barcode;
		.file = x@file;
		.feature = x@feature;
		.peak = x@peak;
		.bmat = x@bmat;
		.pmat = x@pmat;
		.gmat = x@gmat;
		.jmat = x@jmat;
		.smat = x@smat;
		.graph = x@graph;
		.cluster = x@cluster;
		.tsne = x@tsne;
		.umap = x@umap;
		.metaData = x@metaData;
		# a single row or column
       if(!missing(i)){
		   if(max(i) > nrow(x)){
			   stop("idx exceeds number of cells");
		   }
		   if(nrow(.bmat) > 0){.bmat <- .bmat[i,,drop=FALSE]}
		   if(nrow(.pmat) > 0){.pmat <- .pmat[i,,drop=FALSE]}
		   if(nrow(.gmat) > 0){.gmat <- .gmat[i,,drop=FALSE]}	   
		   if(nrow(.jmat@jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
		   if(nrow(.smat@dmat) > 0){.smat <- .smat[i,,drop=FALSE]}
		   if(nrow(.tsne) > 0){.tsne <- .tsne[i,,drop=FALSE]}
		   if(nrow(.umap) > 0){.umap <- .umap[i,,drop=FALSE]}
		   if(nrow(.graph@mat) > 0){.graph <- .graph[i,,drop=FALSE]}
		   if(nrow(.metaData) > 0){.metaData <- .metaData[i,,drop=FALSE]}
		   if(length(.cluster) > 0){.cluster <- .cluster[i,drop=FALSE]}
		   if(length(.barcode) > 0){.barcode <- .barcode[i,drop=FALSE]}
		   if(length(.file) > 0){.file <- .file[i,drop=FALSE]}
	   }
	   if(!missing(j)){
   			mat = match.arg(mat);
	   		if(mat == "bmat"){
	 		   if(ncol(.bmat) > 0){.bmat <- .bmat[,j,drop=FALSE]}
			   if(length(.feature) > 0){.feature <- .feature[j];}	   
	   		}else if(mat == "pmat"){
 	 		   if(ncol(.pmat) > 0){.pmat <- .pmat[,j,drop=FALSE]}
			   if(length(.peak) > 0){.peak <- .peak[j];}	   
	   		}
	   }
	   x@bmat = .bmat;
	   x@pmat = .pmat;
	   x@gmat = .gmat;
	   x@barcode = .barcode;
	   x@file = .file;
	   x@peak = .peak;
	   x@feature = .feature;
	   x@metaData = .metaData;
	   x@umap = .umap;
	   x@feature = .feature;
	   x@jmat = .jmat;
	   x@smat = .smat;
	   x@graph = .graph;
	   x@cluster = .cluster;
	   x@tsne = .tsne;
	   return(x);
})

#' Check barcode existance in snap file
#'
#' This function takes an array of barcodes and a snap-format file as input 
#' and check whether selected barcodes exist in the snap file.
#' 
#' @param barcode An array of selected barcodes.
#' @param file A snap format file.
#' 
#' @return Return an array of logical variable indicates whether the 
#' barcode exists in snap file.
#' 
#' @importFrom rhdf5 h5read
#' @export  
barcodeInSnapFile <- function(barcode, file){
	UseMethod("barcodeInSnapFile", barcode);
}

#' @export  
barcodeInSnapFile.default <- function(barcode, file){
	
	if(missing(file)){
		stop("file is missing");
	}else{
		if(!file.exists(file)){
			stop("file does not exist");
		}
		if(!isSnapFile(file)){
			stop("file is not a snap file");
		}
	}
	
	if(missing(barcode)){
		stop("barcode is missing");
	}else{
		if(class(barcode) != "character"){
			stop("barcode is not a character object");
		}
	}
	
	barcode.ref = as.character(tryCatch(barcode.ref <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
	        options(scipen=999);
	return(barcode %in% barcode.ref);
}


#' Create a snap object from a snap file
#'
#' This function takes a snap-format file as input and create
#' a snap object.
#'
#' @param file A snap-format file name.
#' @param metaData A logical value indicates whether read the meta data [TRUE].
#' @param description Description of the experiment [NULL].
#' @return A snap object
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @export
createSnap <- function(file, metaData, description) {
  UseMethod("createSnap", file);
}

#' @export
createSnap.default <- function(file, metaData=TRUE, description=NULL){	
	
	# close the previously opened H5 file
	H5close();
	# check the input
	if(!file.exists(file)){stop(paste(file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste(file, " is not a snap-format file!", sep=""))};
	if(!(is.logical(metaData))){stop(paste("metaData is not a logical variable!", sep=""))};
		
	
	if(!(is.null(description))){
		if(class(description) != "character"){
			stop("description must be character object")
		}
	}else{
		description=character()
	}
	
	# create an empty snap object
	res = newSnap();
	############################################################################
	message("Epoch: reading the barcode session ...");
	barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @readSnap: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));	
	if(metaData){
		metaData = readMetaData(file);
		if(any((metaData$barcode == barcode) == FALSE)){stop(paste("Error @readSnap: meta data does not match with barcode name!", sep=""))};
		res@metaData = metaData;
	}else{
		metaData = data.frame(barcode=barcode, TN=0, UM=0, PP=0, UQ=0, CM=0);
		res@metaData = metaData;
	}
	nBarcode = length(barcode);
	if((x=nBarcode) == 0L){
		stop("barcode is empty! Does not support reading an empty snap file")
	}
	res@barcode = barcode;
	res@des = description;
	res@file = rep(file, length(res@barcode));
	H5close();
	gc();
	return(res);	
}


#' Add cell-by-bin matrix
#' 
#' This function takes a snap object and snap-format file as input and add the cell-by-bin 
#' matrix to the existing snap object.
#' 
#' @param obj A snap object to add cell-by-bin matrix.
#' @param file A snap file.
#' @param bin.size Cell-by-bin matrix with bin size of bin.size will be added to snap object [5000].
#' @return Return a snap object
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @export
addBmatToSnap <- function(obj, file, bin.size){
  UseMethod("addBmatToSnap", obj);
}

#' @export
addBmatToSnap.default <- function(obj, file, bin.size=5000){	
	# close the previously opened H5 file
	H5close();
	# check the input
	if(class(obj) != "snap"){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))}
	if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
	if(!(bin.size %in% showBinSizes(file))){stop(paste("Error @addBmat: bin.size ", bin.size, " does not exist in ", file, "\n", sep=""))};
	obj@bmat = Matrix(0,0,0);
	message("Epoch: reading cell-bin count matrix session ...");
	############################################################################
	barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));	

	bin.sizeList = showBinSizes(file);
	if(length(bin.sizeList) == 0){stop("Error @addBmat: bin.sizeList is empty! Does not support reading empty snap file")}
	if(!(bin.size %in% bin.sizeList)){stop(paste("Error @addBmat: ", bin.size, " does not exist in bin.sizeList, valid bin.size includes ", toString(bin.sizeList), "\n", sep=""))}
	
	options(scipen=999);
	binChrom = tryCatch(binChrom <- h5read(file, paste("AM", bin.size, "binChrom", sep="/")), error = function(e) {stop(paste("Warning @readaddBmatSnap: 'AM/bin.size/binChrom' not found in ",file))})
	binStart = tryCatch(binStart <- h5read(file, paste("AM", bin.size, "binStart", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})
	if(bin.size == 0){
		binEnd = tryCatch(binEnd <- h5read(file, paste("AM", bin.size, "binEnd", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})		
	}else{
		binEnd = binStart + as.numeric(bin.size) -1
	}
	if((length(binChrom) == 0) || (length(binStart) == 0)){stop("Error @addBmat: bin is empty! Does not support empty snap file")}
	if(length(binChrom) != length(binStart)){
		stop(paste("Error @addBmat: ", "binChrom and binStart has different length!", sep=""))
	}else{
		nBin = length(binChrom);
	}
	bins = GRanges(binChrom, IRanges(as.numeric(binStart),binEnd), name=paste(paste(binChrom, binStart, sep=":"), binEnd, sep="-"));				
	rm(binChrom, binStart);
	obj@feature = bins;
	idx = as.numeric(tryCatch(idx <- h5read(file, paste("AM", bin.size, "idx", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idx' not found in ",file))}));
	idy = as.numeric(tryCatch(idy <- h5read(file, paste("AM", bin.size, "idy", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idy' not found in ",file))}));
	count = as.numeric(tryCatch(count <- h5read(file, paste("AM", bin.size, "count", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/count' not found in ",file))}));	

	if(!all(sapply(list(length(idx),length(idy),length(count)), function(x) x == length(count)))){stop("Error: idx, idy and count has different length in the snap file")}	
	
	ind.sel = which(idx %in% which(barcode %in% obj@barcode));		
	idx = match(idx[ind.sel], sort(unique(idx[ind.sel])));	
	idy = idy[ind.sel];
	count = count[ind.sel];
	nBarcode = length(obj@barcode);
	nBin = length(obj@feature);
	obj@bmat = 	sparseMatrix(i=idx, j =idy, x=count, dims = c(nBarcode, nBin));
	rm(idx, idy, count);
	H5close();
	gc();
	return(obj);
}

#' Add cell-by-peak matrix
#' 
#' This function takes a snap object and snap-format file as input and add the cell-by-peak 
#' matrix to the existing snap object.
#' 
#' @param obj A snap object.
#' @param file A snap file.
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @export
addPmatToSnap <- function(obj, file){
  UseMethod("addPmatToSnap", obj);
}

#' @export
addPmatToSnap.default <- function(obj, file){
        # close the previously opened H5 file
        H5close();
        # check the input
        if(class(obj) != "snap"){stop(paste("Error @addPmat: ", file, " does not exist!", sep=""))}
        if(!file.exists(file)){stop(paste("Error @addPmat: ", file, " does not exist!", sep=""))};
        if(!isSnapFile(file)){stop(paste("Error @addPmat: ", file, " is not a snap-format file!", sep=""))};

        message("Epoch: reading cell-peak count matrix session ...");
        ############################################################################
		barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
        options(scipen=999);
        binChrom = tryCatch(binChrom <- h5read(file, "PM/peakChrom"), error = function(e) {stop(paste("Warning @addPmat: 'PM/peakChrom' not found in ",file))})
        binStart = tryCatch(binStart <- h5read(file, "PM/peakStart"), error = function(e) {stop(paste("Warning @addPmat: 'PM/peakStart' not found in ",file))})
        binEnd = tryCatch(binEnd <- h5read(file, "PM/peakEnd"), error = function(e) {stop(paste("Warning @addPmat: 'PM/peakEnd' not found in ",file))})

        if((length(binChrom) == 0) || (length(binStart) == 0)){stop("Error @readSnap: bin is empty! Does not support empty snap file")}
        if(length(binChrom) != length(binStart)){
                stop(paste("Error @addPmat: ", "binChrom and binStart has different length!", sep=""))
        }else{
                nBin = length(binChrom);
        }

        bins = GRanges(binChrom, IRanges(as.numeric(binStart),binEnd), name=paste(paste(binChrom, binStart, sep=":"), binEnd, sep="-"));
        rm(binChrom, binStart);
        obj@peak = bins;

        idx = as.numeric(tryCatch(idx <- h5read(file, "PM/idx"), error = function(e) {stop(paste("Warning @readSnap: 'PM/idx' not found in ",file))}))
        idy = as.numeric(tryCatch(idy <- h5read(file, "PM/idy"), error = function(e) {stop(paste("Warning @readSnap: 'PM/idy' not found in ",file))}))
        count = as.numeric(tryCatch(count <- h5read(file, "PM/count"), error = function(e) {stop(paste("Warning @readSnap: 'PM/count' not found in ",file))}))

        if(!all(sapply(list(length(idx),length(idy),length(count)), function(x) x == length(count)))){stop("Error: idx, idy and count has different length in the snap file")}
		ind.sel = which(idx %in% which(barcode %in% obj@barcode));		
		idx = match(idx[ind.sel], sort(unique(idx[ind.sel])));	
		idy = idy[ind.sel];
		count = count[ind.sel];
		nBarcode = nrow(obj);
		nPeak = length(obj@peak);
		obj@pmat = 	sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nPeak));
		rm(idx, idy, count);
		H5close();
		gc();
		return(obj);
}

#' Add cell-by-gene matrix
#' 
#' This function takes a snap object and snap-format file as input and add the cell-by-gene 
#' matrix to the existing snap object.
#' 
#' @param obj A snap object to add cell-by-bin matrix.
#' @param file A snap file.
#' @return Return a snap object
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @export
addGmatToSnap <- function(obj, file) {
  UseMethod("addGmatToSnap", obj);
}

#' @export
addGmatToSnap.default <- function(obj, file){
        # close the previously opened H5 file
        H5close();
        # check the input
        if(class(obj) != "snap"){stop(paste("Error @addGmat: ", file, " does not exist!", sep=""))}
        if(!file.exists(file)){stop(paste("Error @addGmat: ", file, " does not exist!", sep=""))};
        if(!isSnapFile(file)){stop(paste("Error @addGmat: ", file, " is not a snap-format file!", sep=""))};

        message("Epoch: reading cell-gene count matrix session ...");
        ############################################################################
		barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
        options(scipen=999);
        geneName = tryCatch(geneName <- h5read(file, "GM/name"), error = function(e) {stop(paste("Warning @addGmat: 'GM/name' not found in ",file))})
        if(length(geneName) == 0){
			stop("Error @addGmat: GM is empty")
		}

        idx = as.numeric(tryCatch(idx <- h5read(file, "GM/idx"), error = function(e) {stop(paste("Warning @addGmat: 'GM/idx' not found in ",file))}))
        idy = as.numeric(tryCatch(idy <- h5read(file, "GM/idy"), error = function(e) {stop(paste("Warning @addGmat: 'GM/idy' not found in ",file))}))
        count = as.numeric(tryCatch(count <- h5read(file, "GM/count"), error = function(e) {stop(paste("Warning @addGmat: 'GM/count' not found in ",file))}))

        if(!all(sapply(list(length(idx),length(idy),length(count)), function(x) x == length(count)))){stop("Error: idx, idy and count has different length in the snap GM session")}		
		ind.sel = which(idx %in% which(barcode %in% obj@barcode));		
		idx = match(idx[ind.sel], sort(unique(idx[ind.sel])));	
		idy = idy[ind.sel];
		count = count[ind.sel];
		nBarcode = nrow(obj);
		nGene = length(geneName);
		obj@gmat = 	sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nGene));
		colnames(obj@gmat) = geneName;
		rm(idx, idy, count);
		H5close();
		gc();
		return(obj);
}

#' Create a snap object from cell-by-bin matrix
#' 
#' This function takes a cell-by-bin count matrix as input and returns a snap object.
#' 
#' @param mat A sparse matrix
#' @param barcodes Corresponding barcodes
#' @param bins A GenomicRanges object for the genomic coordinates of the bins
#' 
#' @return Return a snap object
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @export
createSnapFromBmat <- function(mat, barcodes, bins) {
  UseMethod("createSnapFromBmat", mat);
}

#' @export
createSnapFromBmat.default <- function(mat, barcodes, bins){
	if(missing(mat) || missing(barcodes) || missing(bins)){
		stop("mat or barcodes or bins is missing");
	}

	if(!(class(mat) == "dsCMatrix" || class(mat) == "dgCMatrix")){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(class(bins) != "GRanges"){
		stop("'bins' is not a GRanges object")
	}
	
	if(length(bins) != ncol(mat)){
		stop("'mat' has different number of columns with number of bins");
	}
	
	obj = newSnap();
	obj@bmat = mat;
	obj@barcode = barcodes;
	obj@feature = bins;
	return(obj);
}

#' Create a snap object from cell-by-peak matrix
#' 
#' This function takes a cell-by-peak count matrix as input and returns a snap object.
#' 
#' @param mat A sparse matrix
#' @param barcodes Corresponding barcodes
#' @param peaks A GRanges object for the genomic coordinates of peaks
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @export
createSnapFromPmat <- function(mat, barcodes, peaks) {
  UseMethod("createSnapFromPmat", mat);
}

#' @export
createSnapFromPmat.default <- function(mat, barcodes, peaks){
	if(missing(mat) || missing(barcodes) || missing(peaks)){
		stop("mat or barcodes or peaks is missing");
	}

	if(!(class(mat) == "dsCMatrix" || class(mat) == "dgCMatrix")){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(class(peaks) != "GRanges"){
		stop("'peaks' is not a GRanges object")
	}
	if(length(peaks) != ncol(mat)){
		stop("'mat' has different number of columns with number of peaks");
	}
	
	obj = newSnap();
	obj@pmat = mat;
	obj@barcode = barcodes;
	obj@peak = peaks;
	return(obj);
}

#' Create a snap object from cell-by-gene matrix
#' 
#' This function takes a cell-by-gene count matrix as input and returns a snap object.
#' 
#' @param mat A sparse matrix
#' @param barcodes An array of characters for barcodes
#' @param gene.names An array of characters for gene names
#' @importFrom rhdf5 h5read H5close
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import Matrix
#' @examples
#' gmat = Matrix(sample.int(100, 100), nrow=10, ncol=10, sparse=TRUE);
#' barcodes = paste("barcode", 1:10, sep=".");
#' genes = paste("gene", 1:10, sep=".");
#' createSnapFromGmat(gmat, barcodes, gene.names=genes);
#' @export
createSnapFromGmat <- function(mat, barcodes, gene.names) {
  UseMethod("createSnapFromGmat");
}

#' @export
createSnapFromGmat.default <- function(mat, barcodes, gene.names){
	if(missing(mat) || missing(barcodes) || missing(gene.names)){
		stop("mat or barcodes or gene.names is missing");
	}

	if(!(class(mat) == "dsCMatrix" || class(mat) == "dgCMatrix")){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(class(gene.names) != "character"){
		stop("'gene.names' is not a character object")
	}
	if(length(gene.names) != ncol(mat)){
		stop("'mat' has different number of columns with number of gene.names");
	}
	
	obj = newSnap();
	obj@gmat = mat;
	obj@barcode = barcodes;
	colnames(obj@gmat) = gene.names;
	return(obj);
}

#' Read meta data from a snap file
#'
#' Take a snap file as input and read the barcode session only.
#'
#' @param file character for the snap-format file name which the data are to be read from.
#' @return A data frame contains barcodes and its attributes
#'
#' @importFrom rhdf5 h5read
#' @export
readMetaData <- function(file) {
  UseMethod("readMetaData", file);
}

#' @export
readMetaData.default <- function(file){
	if(!file.exists(file)){stop(paste("Error @readMetaData: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @readMetaData: ", file, " is not a snap-format file!", sep=""))};
	# close the previously opened H5 file
	H5close();
	barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @readSnap: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
	TN = as.numeric(tryCatch(TN <- h5read(file, '/BD/TN'), error = function(e) {print(paste("Warning @readMetaData: 'BD/TN' not found in ",file)); return(c())}));
	UM = as.numeric(tryCatch(UM <- h5read(file, '/BD/UM'), error = function(e) {print(paste("Warning @readMetaData: 'BD/UM' not found in ",file)); return(c())}));
	PP = as.numeric(tryCatch(PP <- h5read(file, '/BD/PP'), error = function(e) {print(paste("Warning @readMetaData: 'BD/PP' not found in ",file)); return(c())}));
	UQ = as.numeric(tryCatch(UQ <- h5read(file, '/BD/UQ'), error = function(e) {print(paste("Warning @readMetaData: 'BD/UQ' not found in ",file)); return(c())}));
	CM = as.numeric(tryCatch(CM <- h5read(file, '/BD/CM'), error = function(e) {print(paste("Warning @readMetaData: 'BD/CM' not found in ",file)); return(c())}));
	data.frame(barcode, TN, UM, PP, UQ, CM)
}


#' Export barcode meta data
#' 
#' This function takes a snap object as input and export its barcode and corresponding attributes
#' 
#' @param obj A snap object
#' @param file Output file name
#' @param slot.names Name of slots to be exported c('barcode', 'tsne', 'umap', 'cluster', 'metaData')
#' @importFrom methods slot
#' @importFrom utils write.table
#' @export
exportMetaData <- function(obj, file, slot.names){
    UseMethod("exportMetaData", obj);	
}

#' @export
exportMetaData.default <- function(obj, file, slot.names=c('barcode', 'cluster', 'tsne', 'umap', 'metaData')){
	subset.names.ref = c('barcode', 'cluster', 'tsne', 'umap', 'metaData');
	
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(class(obj) != "snap"){
			stop("'obj' is not a snap object")
		}
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
	}
	
	if(file.exists(file)){
		stop("file already exists, remove it first")		
	}
	
	if(missing(slot.names)){
		stop("slot.names is missing")
	}

	if(!all(slot.names %in% subset.names.ref)){
		stop("'slot.names' must be subset of c('barcode', 'cluster', 'tsne', 'umap', 'metaData')");		
	}
	
	
	metaData.ls = lapply(as.list(slot.names), function(x){
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


#' Check correlation of cell-by-bin matrix
#'
#' This function takes one or two snap object as input and calculate 
#' the correlation between cell-by-bin matrix between replicates. If
#' obj2 is NULL, obj1 will be split into two pseudo replicates evenly 
#' and the correlaion between these two pseudo-replicates will be 
#' calcualted and returned. For obj1 and obj2, the cell-by-bin matrix
#' must not be empty. This function helps check whether the current 
#' cell-by-bin matrix is sufficient for downstream analysis. If the 
#' pearson correlation is less than 0.95 recommend to use a bigger bin.size. 
#' 
#' @param obj1 A snap object for replicate 1
#' @param obj2 A snap object for replicate 2 [NULL].
#' @return Return pearson correlation between replicates.
#' @importFrom stats cor
#' @export
calBmatCor <- function(obj1, obj2) {
  UseMethod("calBmatCor", obj1);
}

#' @export
calBmatCor.default <- function(obj1, obj2=NULL){
	if(missing(obj1)){
		stop("obj1 is missing.")
	}else{
		if(class(obj1) != "snap"){
			stop("obj1 is not a snap object")
		}
		
		if((x=nrow(obj1@bmat)) == 0){
			stop("cell-by-bin matrix is empty")		
		}
		
		if((x=max(obj1@bmat)) > 1){
			obj1 = makeBinary(obj1, mat="bmat");
		}
	}
	
	if(!is.null(obj2)){
		# check if obj2 is a snap object
		if(class(obj2) != "snap"){
			stop("obj2 is not a snap object")
		}
		
		# check if obj2 has the same features
		if(any(obj1@feature$name != obj2@feature$name)){
			stop("'obj1' and 'obj2' have different features")			
		}
		if(max(obj2@bmat) > 1){
			obj2 = makeBinary(obj2, mat="bmat");									
		}
		cov1 = log(Matrix::colSums(obj1@bmat) + 1, 10);
		cov2 = log(Matrix::colSums(obj2@bmat) + 1, 10);			
	}else{
		ncell = nrow(obj1);
		idx1 = sort(sample(seq(ncell), ncell/2));
		idx2 = setdiff(seq(ncell), idx1);
		cov1 = log(Matrix::colSums(obj1@bmat[idx1,]) + 1, 10);
		cov2 = log(Matrix::colSums(obj1@bmat[idx2,]) + 1, 10);	
	}		
	return(cor(cov1, cov2, method="pearson"));
}


#' Combine snap objects
#'
#' Takes two snap objects and combines them.
#'
#' @param obj1 a snap object
#' @param obj2 a snap object
#' @return a combined snap object
#' @export
snapRbind <- function(obj1, obj2){
  UseMethod("snapRbind", obj1);
}

#' @export
snapRbind.default <- function(obj1, obj2){
	if(!is.snap(obj1)){stop(paste("Error @snapRbind: obj1 is not a snap object!", sep=""))};
	if(!is.snap(obj2)){stop(paste("Error @snapRbind: obj2 is not a snap object!", sep=""))};

	# barcode from obj1 and obj2
	barcode1 = obj1@barcode;
	barcode2 = obj2@barcode;	
	
	# check barcode name, if there exists duplicate barcode raise error and exist
	if(length(unique(c(barcode1, barcode2))) < length(barcode1) + length(barcode2)){
		stop("Error: @snapRbind: identifcal barcodes found in obj1 and obj2!")
	}
	barcode = c(barcode1, barcode2);
	
	# check meta data
	if(nrow(obj1@metaData) > 0 && nrow(obj2@metaData) > 0){
		metaData = rbind(obj1@metaData, obj2@metaData);		
	}else{
		metaData = data.frame();
	}
	
	# check feature
	feature1 = obj1@feature;
	feature2 = obj2@feature;
	if((length(feature1) == 0) != (length(feature2) == 0)){
		stop("different feature found in obj1 and obj2!")
	}else{
		if(length(feature1) > 0){
			if(FALSE %in% (feature1$name == feature2$name)){
				stop("Error: @snapRbind: different feature found in obj1 and obj2!")
			}
			feature = feature1;					
		}else{
			feature = feature1;								
		}
	}
	
	# check peak
	peak1 = obj1@peak;
	peak2 = obj2@peak;
	if((length(peak1) == 0) != (length(peak2) == 0)){
		stop("different peak found in obj1 and obj2!")
	}else{
		if(length(peak1) > 0){
			if(FALSE %in% (peak1$name == peak2$name)){
				stop("Error: @snapRbind: different feature found in obj1 and obj2!")
			}
			peak = peak1;					
		}else{
			peak = peak1;								
		}
	}
	
	# check bmat	
	bmat1 = obj1@bmat;
	bmat2 = obj2@bmat;
	if((length(bmat1) == 0) != (length(bmat2) == 0)){
		stop("bmat has different dimentions in obj1 and obj2!")
	}else{
		bmat = Matrix::rBind(bmat1, bmat2);
	}

	# check gmat	
	gmat1 = obj1@gmat;
	gmat2 = obj2@gmat;
	if((length(gmat1) == 0) != (length(gmat2) == 0)){
		stop("gmat has different dimentions in obj1 and obj2!")
	}else{
		gmat = Matrix::rBind(gmat1, gmat2);
	}

	# check pmat	
	pmat1 = obj1@pmat;
	pmat2 = obj2@pmat;
	if((length(pmat1) == 0) != (length(pmat2) == 0)){
		stop("pmat has different dimentions in obj1 and obj2!")
	}else{
		pmat = Matrix::rBind(pmat1, pmat2);
	}
	gc();
	res = newSnap();
	res@barcode = barcode;
	res@metaData = metaData;
	res@bmat = bmat;
	res@pmat = pmat;
	res@feature = feature;
	res@peak = peak;
	res@gmat = gmat;
	return(res)
}


#' Cell filtration
#'
#' This function takes a snap object as input and filter cells based on given cutoffs. 
#' We next identify the high-quality barcode based on the following metrices: 
#' 
#' 1) fragment.num - total number of fragments per barcode;
#' 2) UMI - unique molecular identifier;
#' 3) mito.ratio - mitochondrial ratio;
#' 4) dup.ratio - PCR duplicate ratio;
#' 5) pair.ratio - properly paired ratio;
#' 6) umap.ratio - uniquely mapped ratio;
#' 
#' Note we no longer use reads in peak ratio as a metric for cell selection mainly for two reasons. 
#' Reads-in-peak ration is highly cell type specific. For instance, according to published single 
#' cell ATAC-seq, human fibroblast (BJ) cells have significantly higher reads in peak ratio (40-60%) 
#' versus 20-40% for GM12878 cells. Similarly, in mammalian brain, glia cells overall have very 
#' different reads in peak ratio distribution compared to neuronal cells. We suspect this may reflect 
#' the nucleus size or global chromatin accessibility. Second, pre-defined set of accessibility peaks 
#' are incomplete and biased to the dominant populations. 
#' 
#' @param obj A snap object.
#' @param subset.names Attributes used to filter cells c('fragment.num', 'UMI', 'mito.ratio', 'umap.ratio', 'dup.ratio', 'pair.ratio'). 
#' @param low.thresholds Low cutoffs for the parameters (default is -Inf)
#' @param high.thresholds High cutoffs for the parameters (default is Inf)
#'
#' @return Returns a snap object containing only the relevant subset of cells
#' 
#' @export
filterCells <- function(obj, subset.names, low.thresholds, high.thresholds) {
  UseMethod("filterCells", obj);
}

#' @export
filterCells.default <- function(
	obj, 
	subset.names, 
	low.thresholds, 
	high.thresholds
){
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(class(obj) != "snap"){
			stop("obj is not a snap object");
		}
		metaData = obj@metaData;		
		if((x=nrow(metaData))==0L){
			stop("metaData is empty")
		}
	}
	
    if (missing(x = low.thresholds)){
      low.thresholds <- replicate(n = length(x = subset.names), expr = -Inf);
    }

    if (missing(x = high.thresholds)) {
      high.thresholds <- replicate(n = length(x = subset.names), expr = Inf);
    }

	if(!all(subset.names %in% c('fragment.num', 'mito.ratio', 'umap.ratio', 'dup.ratio', 'pair.ratio', 'UMI'))){
		stop("'subset.names' must be subset of c('UMI', 'fragment.num', 'mito.ratio', 'umap.ratio', 'dup.ratio', 'pair.ratio')");		
	}

    length.check <- sapply(
      X = list(subset.names, low.thresholds, high.thresholds),
      FUN = length
    )
	
    if (length(x = unique(x = length.check)) != 1L) {
      stop("'subset.names', 'low.thresholds', and 'high.thresholds' must all have the same length");
    }

	data.subsets <- data.frame(name=subset.names, lw=low.thresholds, hg=high.thresholds);
	
	metaData$UMI = metaData$UQ;
	metaData$fragment.num = metaData$TN;
	metaData$mito.ratio = (metaData$CM)/(metaData$UQ);
	metaData$umap.ratio = (metaData$UM)/(metaData$TN);
	metaData$dup.ratio = 1 - (metaData$UQ)/(metaData$PP);
	metaData$pair.ratio = (metaData$PP)/(metaData$UM);
	for(i in seq(nrow(data.subsets))){		
		f = as.character(data.subsets$name[i]);
		names.use <- which(colnames(metaData) %in% f);
		idx = which((metaData[,names.use] >= data.subsets$lw[i]) & (metaData[,names.use] <= data.subsets$hg[i]));
		obj = obj[idx,];
		metaData = metaData[idx,]
	}
	return(obj);
}

#' Check a snap-format file
#' 
#' This function takes a file name as input and check if the file is a snap-formated file 
#' @param file A file name.
#' @return Return a logical variable indicates whether file is a snap file.
#' @importFrom rhdf5 h5read h5ls
#' 
#' @export
#' 
isSnapFile <- function(file) {
  UseMethod("isSnapFile", file);
}

#' @export
isSnapFile.default <- function(file){
	if(!file.exists(file)){ 
		stop("file does not exist!")
	}

	monals = tryCatch(
		monals <- h5ls(file), 
		error = function(e) {
			return(character())
	})
	
	if(length(monals) != 0){
		magicString = as.character(tryCatch(magicString <- h5read(file, '/HD/MG'), error = function(e) {print(paste("Warning @isSnapFile: 'HD/MG' not found in ", file)); return(character())}))
		if(length(magicString)==0){
			return(FALSE);
		}
		if(magicString == "SNAP" || magicString == "MONA"){
			return(TRUE);
		}
	}
	return(FALSE);
}

#' Calculate cell-by-gene count matrix
#'
#' This function takes a snap obj and converts the cell-by-bin 
#' matrix into the cell-by-gene count matrix.
#'
#' @param obj A snap object.
#' @param gene A GRanges object that contains the genomic intervals of genes.
#' @param num.cores Number of CPU processers for computation [1].
#' @param mat Matrix slot is used to calculatd the cell-by-gene matrix c("bmat", "pmat").
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom parallel mclapply
#' @importFrom methods slot
#' @return Returns a Snap obj with the cell-by-gene matrix stored in obj@gmat
#' 
#' @export
calCellGeneTable <- function(obj, gene, num.cores, mat){
  UseMethod("calCellGeneTable", obj);
}

#' @export
calCellGeneTable.default <- function(
	obj, 
	gene, 
	num.cores=1, 
	mat=c("bmat", "pmat")
){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(class(obj) != "snap"){
			stop("obj is not a snap obj")
		}
	}
	mat = match.arg(mat);
	data.use = methods::slot(obj, mat);
	if((x=nrow(data.use)) == 0L){
		stop("mat is empty")
	}
	
	if(mat == "bmat"){
		feature.use = obj@feature;
	}else if(mat == "pmat"){
		feature.use = obj@peak;
	}
	
	if((x=length(feature.use)) == 0L){
		stop("column feature for mat is empty")
	}
	
	if(missing(gene)){
		stop("gene is missing");
	}else{
		if(class(gene) != "GRanges"){
			stop("'gene' is not a GenomicRanges obj");
		}
		if((x==length(gene)) == 0L){
			stop("gene is empty");
		}		
	}
		
	if(any(duplicated(gene$gene_name))){
		stop("'gene' contains duplicate gene names, remove duplicate first");
	}

	# find overlap between feature and genes
	ov = data.frame(findOverlaps(feature.use, gene));
	if(nrow(ov) > 0){
		# calculate gene count vector per cell in parallel;
		ov.ls <- split(ov, ov$subjectHits);
		# generate the count vector per cell in parallel;
		count.ls <- mclapply(ov.ls, function(x){
			if(nrow(x) > 1){
				return(Matrix::rowSums(mat[,x[,1]]))
			}else{
				return(mat[,x[,1]])
			}
		}, mc.cores=num.cores);
		# combine vectors and create a count matrix;
		count.mt = t(as.matrix(do.call(rbind, count.ls)));
		count_table = Matrix(0, nrow=nrow(obj), ncol=length(gene), sparse=TRUE);
		count_table[,as.numeric(names(ov.ls))] = count.mt;
		count_table = count_table;		
	}else{
		count_table = Matrix(0, nrow(obj), ncol=length(gene), sparse=TRUE);
	}
	obj@gmat = count_table;
	colnames(obj@gmat) = gene$gene_name;
	return(obj);
}

