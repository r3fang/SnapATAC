#' Create an empty snap object
#'
#' This function creates an empty snap object
#'
#' @examples
#' x.sp = newSnap();
#' 
#' @return An empty snap object
#'
#' @importFrom methods new
#' @importFrom GenomicRanges GRanges
#' @export
#' 
newSnap <- function () {
	metaData=data.frame();
	des = character()
	file = as.character(c());
	sample = as.character(c());
	barcode = as.character(c());
	feature = GRanges();
	peak = GRanges();		
	metaData = data.frame();
	bmat=Matrix(nrow=0, ncol=0, sparse=TRUE);		
	pmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	gmat=Matrix(nrow=0, ncol=0, sparse=TRUE);
	mmat=matrix(0,0,0);
	jmat=newJaccard();
	smat=newDimReduct();	
	graph=newKgraph();
	regModel=c();
	tsne=matrix(nrow=0, ncol=0);	
	umap=matrix(nrow=0, ncol=0);	
	cluster=factor();
	res = new("snap", 
			  des=des,
			  file=file,
			  sample=sample,
			  barcode=barcode, 
			  feature=feature, 
			  peak=peak, 
			  metaData=metaData, 
			  bmat=bmat, 
			  pmat=pmat, 
			  gmat=gmat,
			  mmat=mmat, 
			  jmat=jmat, 
			  smat=smat, 
			  graph=graph, 
			  tsne=tsne, 
			  umap=umap, 
			  cluster=cluster
			  );	
}

#' Show bin sizes in a snap file
#'
#' This function takes a snap-format file name as input and check the bin 
#' sizes or resolutions have been generated for count matrix
#'
#' @param file character. input snap-format file name
#'
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' showBinSizes(file.name);
#' 
#' @return integer vector. A vector of integers indicating the bin sizes
#'
#' @export
#'
showBinSizes <- function(file) {
  UseMethod("showBinSizes", file);
}

#' @export
showBinSizes.default <- function(file){
	# this is to check what are the binSizes have been generated
	if(!file.exists(file)){stop("file does not exist!")}
	if(!isSnapFile(file)){stop("input file is not a snap file\n")}		
	binSizeList = tryCatch(binSizeList <- h5read(file, '/AM/binSizeList'), error = function(e) {print(paste("Warning @check.bin.size: '/AM/binSizeList' not found in ", file)); return(numeric())})
	return(binSizeList);
}

#' summarySnap for snap object.
#'
#' This function takes a snap object and returns the summary statistics
#' @name summarySnap
#' @param obj snap; a snap object
#' @rdname summarySnap-methods
#' @importFrom stats median
#' @exportMethod summarySnap
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' x.sp = createSnap(file.name, sample="demo");
#' summarySnap(x.sp);
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
#' @examples
#' data(demo.sp);
#' nrow(demo.sp);
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
#' @examples
#' data(demo.sp);
#' a = colSums(demo.sp, mat="bmat");
#' b = colSums(demo.sp, mat="pmat");
#' d = colSums(demo.sp, mat="gmat");
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
#' 
#' @examples
#' data(demo.sp);
#' a = rowSums(demo.sp, mat="bmat");
#' b = rowSums(demo.sp, mat="pmat");
#' d = rowSums(demo.sp, mat="gmat");
#' 
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
#'
#' @examples
#' data(demo.sp);
#' a = rowMeans(demo.sp, mat="bmat");
#' b = rowMeans(demo.sp, mat="pmat");
#' d = rowMeans(demo.sp, mat="gmat");
#' 
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
#' 
#' @examples
#' data(demo.sp);
#' a = colMeans(demo.sp, mat="bmat");
#' b = colMeans(demo.sp, mat="pmat");
#' d = colMeans(demo.sp, mat="gmat");
#' 
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
#' @examples
#' data(demo.sp);
#' is.snap(demo.sp);
#' @rdname is.snap-methods
#' @exportMethod is.snap
setGeneric("is.snap", function(obj) standardGeneric("is.snap"))

#' @rdname is.snap-methods
#' @aliases is.snap,snap-method
setMethod("is.snap", "snap", function(obj) return(is(obj, "snap")));

#' subsetting for snap objects
#'
#' This function takes a snap object and returns the subset of snap object.
#' @param x snap; a snap object
#' @param i integer; selected barcode index
#' @param j integer; selected feature index
#' @param mat character; indicates the slot to subsetting
#' @param drop character; 
#' @examples
#' data(demo.sp);
#' demo.sp[1:10,];
#' demo.sp[,1:10,mat="bmat"];
#' demo.sp[,1:10,mat="pmat"];
#' @export
setMethod("[", "snap",
	function(x,i,j,mat=c("bmat", "pmat", "gmat"), drop="missing"){
		.barcode = x@barcode;
		.file = x@file;
		.sample = x@sample;
		.feature = x@feature;
		.peak = x@peak;
		.bmat = x@bmat;
		.pmat = x@pmat;
		.gmat = x@gmat;
		.mmat = x@mmat;
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
		   if(nrow(.mmat) > 0){.mmat <- .mmat[i,,drop=FALSE]}	   
		   if(nrow(.jmat@jmat) > 0){.jmat <- .jmat[i,,drop=FALSE]}
		   if(nrow(.smat@dmat) > 0){.smat <- .smat[i,,drop=FALSE]}
		   if(nrow(.tsne) > 0){.tsne <- .tsne[i,,drop=FALSE]}
		   if(nrow(.umap) > 0){.umap <- .umap[i,,drop=FALSE]}
		   if(nrow(.graph@mat) > 0){.graph <- .graph[i,,drop=FALSE]}
		   if(nrow(.metaData) > 0){.metaData <- .metaData[i,,drop=FALSE]}
		   if(length(.cluster) > 0){.cluster <- .cluster[i,drop=FALSE]}
		   if(length(.barcode) > 0){.barcode <- .barcode[i,drop=FALSE]}
		   if(length(.file) > 0){.file <- .file[i,drop=FALSE]}
		   if(length(.sample) > 0){.sample <- .sample[i,drop=FALSE]}
	   }
	   if(!missing(j)){
   			mat = match.arg(mat);
	   		if(mat == "bmat"){
	 		   if(ncol(.bmat) > 0){.bmat <- .bmat[,j,drop=FALSE]}
			   if(length(.feature) > 0){.feature <- .feature[j];}	   
	   		}else if(mat == "pmat"){
 	 		   if(ncol(.pmat) > 0){.pmat <- .pmat[,j,drop=FALSE]}
			   if(length(.peak) > 0){.peak <- .peak[j];}	   
	   		}else if(mat == "gmat"){
 	 		   if(ncol(.gmat) > 0){.gmat <- .gmat[,j,drop=FALSE]}
	   		}
	   }
	   x@bmat = .bmat;
	   x@pmat = .pmat;
	   x@gmat = .gmat;
	   x@mmat = .mmat;
	   x@barcode = .barcode;
	   x@file = .file;
	   x@sample = .sample;
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
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' barcodes = c("ACATTGGCAACCAGGTTGCTGGTATTGGAAGT", "ACATTGGCAAGAGGCAACAAGGATATCTGAGT");
#' barcodeInSnapFile(barcodes, file.name);
#' 
#' @return Return an array of logical variable indicates whether the 
#' barcode exists in snap file.
#' @importFrom rhdf5 h5read
#' @importFrom methods is
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
		if(!is(barcode, "character")){
			stop("barcode is not a character object");
		}
	}
	
	barcode.ref = as.character(tryCatch(barcode.ref <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @barcodeInSnapFile: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
	        options(scipen=999);
	return(barcode %in% barcode.ref);

}

#' Create a snap object from a snap file
#'
#' This function takes a snap-format file as input and create
#' a snap object.
#'
#' @param file Name of a snap-format file.
#' @param sample A short sample name (i.g. "MOS.rep1").
#' @param description Description of the experiment [NULL].
#' @param do.par A logical variable indicates if run this using multiple processors [FALSE].
#' @param num.cores Number of processers to use [1].
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' demo.sp = createSnap(file.name, sample="demo", do.par=FALSE);
#' 
#' @return A snap object
#' @importFrom rhdf5 h5read
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
createSnap <- function(file, sample, description, do.par, num.cores) {
  UseMethod("createSnap", file);
}

#' @export
createSnap.default <- function(file, sample, description=NULL, do.par=FALSE, num.cores=1){
	
	if(missing(file)){
		stop("file is missing");
	}else{
		if(!is(file, "character")){
			stop("file is not character")
		}
		if(any(duplicated(file))){
			stop("file has duplicate name");
		}		
	}
	
	if(!is.numeric(num.cores)){
		stop("num.cores is not an integer")
	}		
	num.cores = round(num.cores);

	fileList = as.list(file);	
	if(missing(sample)){
		stop("sample is missing");
	}else{
		if(!is(sample, "character")){
			stop("sample is not character")
		}
		if(any(duplicated(sample))){
			stop("sample has duplicate name");
		}
	}
	
	if(length(sample) != length(file)){
		stop("sample has different length with file");
	}
	
	sampleList = as.list(sample);
	
	if(!(is.null(description))){
		if(!is(description, "character")){
			stop("description must be character object")
		}
	}else{
		description=character()
	}
	
	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	message("Epoch: reading the barcode session ...");
	if(do.par){
		obj.ls = mclapply(as.list(seq(fileList)), function(i){
			createSnapSingle(file=fileList[[i]], sample=sampleList[[i]]);
		}, mc.cores=num.cores);
	}else{
		obj.ls = lapply(as.list(seq(fileList)), function(i){
			createSnapSingle(file=fileList[[i]], sample=sampleList[[i]]);
		});		
	}
	
	obj = Reduce(snapRbind, obj.ls);
	rm(obj.ls);
	gc();
	obj@des = description;
	return(obj);
}

#' Add cell-by-bin matrix
#' 
#' This function takes a snap object as input and add the cell-by-bin 
#' matrix to the existing snap object.
#' 
#' @param obj A snap object to add cell-by-bin matrix.
#' @param bin.size Cell-by-bin matrix with bin size of bin.size 
#' will be added to snap object.
#' @param do.par A logical variable indicates whether use multiple processors [FALSE].
#' @param num.cores Number of processors to use [1].
#' 
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' demo.sp = createSnap(file.name, sample="demo", do.par=FALSE);
#' showBinSizes(file.name);
#' demo.sp = addBmatToSnap(demo.sp, bin.size=100000, do.par=FALSE);
#' 
#' @return Return a snap object
#' @importFrom rhdf5 h5read
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
addBmatToSnap <- function(obj, bin.size, do.par, num.cores){
  UseMethod("addBmatToSnap", obj);
}

#' @export
addBmatToSnap.default <- function(obj, bin.size=5000, do.par=FALSE, num.cores=1){	
	# close the previously opened H5 file
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
	}
	
	if(!is.numeric(num.cores)){
		stop("num.cores is not an integer")
	}		
	num.cores = round(num.cores);
	
	fileList = as.list(unique(obj@file));

	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if BM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "AM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "AM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following nsap files do not contain AM session")
		print(fileList[idx])
		stop()
	}
	
	if(any(do.call(c, lapply(fileList, function(x){(bin.size %in% showBinSizes(x))})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){(bin.size %in% showBinSizes(x))})) == FALSE)
		print("error: chosen bin size does not exist in the following snap files")
		print(fileList[idx])
		stop()
	}

	# check if bins match
	bin.list = lapply(fileList, function(x){
		readBins(x, bin.size=bin.size)
	})
	
	if(!all(sapply(bin.list, FUN = identical, bin.list[[1]]))){
		stop("bins does not match between snap files, please regenerate the cell-by-bin matrix by snaptools")
	}
	
	# read the snap object
	message("Epoch: reading cell-bin count matrix session ...");
	if(do.par){
		obj.ls = mclapply(fileList, function(file){
			idx = which(obj@file == file)
			addBmatToSnapSingle(obj[idx,], file, bin.size=bin.size);
		}, mc.cores=num.cores);		
	}else{
		obj.ls = lapply(fileList, function(file){
			idx = which(obj@file == file)
			addBmatToSnapSingle(obj[idx,], file, bin.size=bin.size);
		});		
	}
	
	# combine
	if((x=length(obj.ls)) == 1L){
		res = obj.ls[[1]]
	}else{
		res = Reduce(snapRbind, obj.ls);		
	}
	obj@feature = res@feature;
	obj@bmat = res@bmat;
	rm(res, obj.ls);
	gc()
	return(obj);
}

#' Add cell-by-peak matrix
#' 
#' This function takes a snap object as input and add the 
#' cell-by-peak matrix to the existing snap object.
#' 
#' @param obj A snap object.
#' @param do.par A logical varaible indicates whether to use multiple processors [FALSE].
#' @param num.cores Number of processors to use [1].
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' demo.sp = createSnap(file.name, sample="demo", do.par=FALSE);
#' showBinSizes(file.name);
#' demo.sp = addPmatToSnap(demo.sp, do.par=FALSE);
#' 
#' @importFrom rhdf5 h5read
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @export
addPmatToSnap <- function(obj, do.par, num.cores){
  UseMethod("addPmatToSnap", obj);
}

#' @export
addPmatToSnap.default <- function(obj, do.par=FALSE, num.cores=1){
	# close the previously opened H5 file
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}

	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
	}
	
	if(!is.numeric(num.cores)){
		stop("num.cores is not an integer")
	}		
	num.cores = round(num.cores);
	
	
	fileList = as.list(unique(obj@file));

	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if PM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "PM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "PM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following snap files do not contain PM session")
		print(fileList[idx])
		stop()
	}
	
	# check if bins match
	peak.list = lapply(fileList, function(x){
		readPeaks(x)
	})
	
	if(!all(sapply(peak.list, FUN = identical, peak.list[[1]]))){
		stop("peaks does not match between snap files, please regenerate the cell-by-peak matrix by snaptools using the same peak file")
	}
	
	# read the snap object
	message("Epoch: reading cell-peak count matrix session ...");
	if(do.par){
		obj.ls = mclapply(fileList, function(file){
			idx = which(obj@file == file)
			addPmatToSnapSingle(obj[idx,], file);
		}, mc.cores=num.cores);		
	}else{
		obj.ls = lapply(fileList, function(file){
			idx = which(obj@file == file)
			addPmatToSnapSingle(obj[idx,], file);
		});				
	}
	
	# combine
	if((x=length(obj.ls)) == 1L){
		res = obj.ls[[1]]
	}else{
		res = Reduce(snapRbind, obj.ls);		
	}
	
	# re-order the matrix
	o1 = paste(obj@file, obj@barcode, sep=".");
	o2 = paste(res@file, res@barcode, sep=".");
	obj@peak = res@peak;
	obj@pmat = res@pmat[match(o1, o2),];
	rm(obj.ls, res, o1, o2);
	gc();
	return(obj);
}

#' Add cell-by-gene matrix
#' 
#' This function takes a snap object as input and add the cell-by-gene 
#' matrix to the existing snap object.
#' 
#' @param obj A snap object to add cell-by-bin matrix.
#' @param do.par A logical variable indicates whether to use multiple processors [FALSE].
#' @param num.cores Number of processors to use.
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' demo.sp = createSnap(file.name, sample="demo", do.par=FALSE);
#' demo.sp = addGmatToSnap(demo.sp, do.par=FALSE);
#' 
#' @return Return a snap object
#' @importFrom rhdf5 h5read 
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @export
addGmatToSnap <- function(obj, do.par, num.cores) {
  UseMethod("addGmatToSnap", obj);
}

#' @export
addGmatToSnap.default <- function(obj, do.par=FALSE, num.cores=1){	
	# close the previously opened H5 file
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
	}
	
	if(!is.numeric(num.cores)){
		stop("num.cores is not an integer")
	}		
	num.cores = round(num.cores);
	
	fileList = as.list(unique(obj@file));

	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if GM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "GM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "GM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following nsap files do not contain GM session")
		print(fileList[idx])
		stop()
	}
	
	# read the snap object
	message("Epoch: reading cell-gene count matrix session ...");
	if(do.par){
		obj.ls = mclapply(fileList, function(file){
			idx = which(obj@file == file);
			addGmatToSnapSingle(obj[idx,], file);
		}, mc.cores=num.cores);		
	}else{
		obj.ls = lapply(fileList, function(file){
			idx = which(obj@file == file);
			addGmatToSnapSingle(obj[idx,], file);
		});				
	}

	# combine
	if((x=length(obj.ls)) == 1L){
		res = obj.ls[[1]]
	}else{
		res = Reduce(snapRbind, obj.ls);		
	}
	o1 = paste(obj@file, obj@barcode, sep=".");
	o2 = paste(res@file, res@barcode, sep=".");
	obj@gmat = res@gmat[match(o1, o2),];
	rm(obj.ls, res, o1, o2);
	gc();
	return(obj);
}

#' Remove cell-by-bin matrix
#' 
#' This function takes a snap object as input and removes the cell-by-bin 
#' matrix in the existing snap object. Report error when cell-by-bin matrix is empty.
#' 
#' @param obj A snap object to remove cell-by-bin matrix.
#' @examples
#' data(demo.sp)
#' rmBmatFromSnap(demo.sp)
#' 
#' @return Return a snap object without cell-by-bin matrix 
#' @importFrom methods slot
#' @importFrom GenomicRanges GRanges
#' @import Matrix
#' @export
rmBmatFromSnap <- function(obj){
  UseMethod("rmBmatFromSnap", obj);
}

#' @export
rmBmatFromSnap.default <- function(obj){	
	# close the previously opened H5 file
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
		data.use = methods::slot(obj, "bmat");
		if((x=nrow(data.use)) == 0L){
			stop("cell-by-bin matrix does not exist in obj");
		}
	}
	
	obj@bmat = Matrix(0,0,0, sparse=TRUE);
	obj@feature = GRanges();
	return(obj);
}


#' Remove cell-by-peak matrix
#' 
#' This function takes a snap object as input and removes the cell-by-peak 
#' matrix in the existing snap object. Report error when cell-by-peak matrix is empty.
#' 
#' @param obj A snap object to remove cell-by-peak matrix.
#' @examples
#' data(demo.sp)
#' rmPmatFromSnap(demo.sp)
#' 
#' @return Return a snap object without cell-by-peak matrix
#' @importFrom methods slot
#' @importFrom GenomicRanges GRanges
#' @import Matrix
#' @export
rmPmatFromSnap <- function(obj){
  UseMethod("rmPmatFromSnap", obj);
}

#' @export
rmPmatFromSnap.default <- function(obj){	
	# close the previously opened H5 file
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
		data.use = methods::slot(obj, "pmat");
		if((x=nrow(data.use)) == 0L){
			stop("cell-by-peak matrix does not exist in obj");
		}
	}
	
	obj@pmat = Matrix(0,0,0, sparse=TRUE);
	obj@peak = GRanges();
	return(obj);
}

#' Remove cell-by-gene matrix
#' 
#' This function takes a snap object as input and removes the cell-by-gene 
#' matrix in the existing snap object. Report error when cell-by-gene matrix 
#' is empty.
#' 
#' @param obj A snap object to remove cell-by-gene matrix.
#' @examples
#' data(demo.sp)
#' rmPmatFromSnap(demo.sp)
#' 
#' @return Return a snap object without cell-by-peak matrix
#' @importFrom methods slot
#' @importFrom GenomicRanges GRanges
#' @import Matrix
#' @export
rmGmatFromSnap <- function(obj){
  UseMethod("rmGmatFromSnap", obj);
}

#' @export
rmGmatFromSnap.default <- function(obj){	
	# close the previously opened H5 file
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
		data.use = methods::slot(obj, "gmat");
		if((x=nrow(data.use)) == 0L){
			stop("cell-by-gene matrix does not exist in obj");
		}
	}
	
	obj@gmat = Matrix(0,0,0, sparse=TRUE);
	return(obj);
}

#' Create a snap object from cell-by-bin matrix
#' 
#' This function takes a cell-by-bin count matrix as input and returns a snap object.
#' 
#' @param mat A sparse matrix
#' @param barcodes Corresponding barcodes
#' @param bins A GenomicRanges object for the genomic coordinates of the bins
#' @examples
#' library("GenomicRanges");
#' mat = Matrix(sample(0:10, 100, replace=TRUE),sparse=TRUE, ncol=5);
#' barcodes = paste("barcode", seq(nrow(mat)), sep=".");
#' chroms = c("chr1", "chr1", "chr1", "chr1", "chr1");
#' pos = c(1, 5001, 10001, 15001, 20001);
#' bins = GRanges(chroms, IRanges(pos, pos+5000));
#' x.sp = createSnapFromBmat(
#'	mat, 
#'	barcodes=barcodes,
#'	bins=bins
#'	);
#' @return Return a snap object
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @export
createSnapFromBmat <- function(mat, barcodes, bins) {
  UseMethod("createSnapFromBmat", mat);
}

#' @export
createSnapFromBmat.default <- function(mat, barcodes, bins){
	if(missing(mat) || missing(barcodes) || missing(bins)){
		stop("mat or barcodes or bins is missing");
	}

	if(!(is(mat, "dsCMatrix") || is(mat, "dgCMatrix"))){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(!is(bins, "GRanges")){
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
#' @examples
#' library("GenomicRanges");
#' mat = Matrix(sample(0:10, 100, replace=TRUE),sparse=TRUE, ncol=5);
#' barcodes = paste("barcode", seq(nrow(mat)), sep=".");
#' chroms = c("chr1", "chr1", "chr1", "chr1", "chr1");
#' pos = c(1, 5001, 10001, 15001, 20001);
#' peaks = GRanges(chroms, IRanges(pos, pos+100));
#' x.sp = createSnapFromPmat(
#'	mat, 
#'	barcodes=barcodes,
#'	peaks=peaks
#'	);
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @export
createSnapFromPmat <- function(mat, barcodes, peaks) {
  UseMethod("createSnapFromPmat", mat);
}

#' @export
createSnapFromPmat.default <- function(mat, barcodes, peaks){
	if(missing(mat) || missing(barcodes) || missing(peaks)){
		stop("mat or barcodes or peaks is missing");
	}

	if(!(is(mat, "dsCMatrix") || is(mat, "dgCMatrix"))){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(!is(peaks, "GRanges")){
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
#' @examples
#' mat = Matrix(sample(0:10, 100, replace=TRUE),sparse=TRUE, ncol=5);
#' barcodes = paste("barcode", seq(nrow(mat)), sep=".");
#' gene.names = paste("genes", seq(ncol(mat)), sep=".");
#' x.sp = createSnapFromGmat(
#'	mat, 
#'	barcodes=barcodes,
#'	gene.names=gene.names
#'	);
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @export
createSnapFromGmat <- function(mat, barcodes, gene.names) {
  UseMethod("createSnapFromGmat");
}

#' @export
createSnapFromGmat.default <- function(mat, barcodes, gene.names){
	if(missing(mat) || missing(barcodes) || missing(gene.names)){
		stop("mat or barcodes or gene.names is missing");
	}

	if(!(is(mat, "dsCMatrix") || is(mat, "dgCMatrix"))){
		stop("'mat' is not a sparse matrix");
	}

	if(length(barcodes) != nrow(mat)){
		stop("'mat' has different number of rows with number of barcodes");
	}
	
	if(!is(gene.names, "character")){
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
#'
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' md = readMetaData(file.name);
#' 
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
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
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
#' @importFrom methods is
#' @examples
#' data(demo.sp);
#' exportMetaData(demo.sp, file="demo.metadata.txt", slot.names=c("barcode", "tsne"));
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
		if(!is(obj, "snap")){
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
#' obj2 is NULL, obj1 will be randomly split into two pseudo replicates  
#' and the correlaion between these two pseudo-replicates will be 
#' calcualted and returned. For obj1, the cell-by-bin matrix
#' cannot be empty. This function helps check whether the current 
#' cell-by-bin matrix is sufficient for downstream analysis. If the 
#' pearson correlation is less than 0.95 recommend to use a bigger bin.size. 
#' 
#' @param obj1 A snap object for replicate 1
#' @param obj2 A snap object for replicate 2 [NULL].
#' @examples
#' data(demo.sp);
#' calBmatCor(demo.sp)
#' 
#' @return Return pearson correlation between replicates.
#' @importFrom stats cor
#' @importFrom methods is
#' @export
calBmatCor <- function(obj1, obj2) {
  UseMethod("calBmatCor", obj1);
}

#' @export
calBmatCor.default <- function(obj1, obj2=NULL){
	if(missing(obj1)){
		stop("obj1 is missing.")
	}else{
		if(!is(obj1, "snap")){
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
		if(!is(obj2, "snap")){
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
	if(!is.snap(obj1)){stop(paste("Error @snapRbind: obj1 is not a snap object!", sep=""))};
	if(!is.snap(obj2)){stop(paste("Error @snapRbind: obj2 is not a snap object!", sep=""))};

	# barcode from obj1 and obj2
	barcode1 = paste(obj1@file, obj1@barcode, sep=".");
	barcode2 = paste(obj2@file, obj2@barcode, sep=".");	
	
	# check barcode name, if there exists duplicate barcode raise error and exist
	if(length(unique(c(barcode1, barcode2))) < length(barcode1) + length(barcode2)){
		stop("Error: @snapRbind: identifcal barcodes found in obj1 and obj2!")
	}
	rm(barcode1, barcode2)
	gc()
	
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
	gc()
	
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
	rm(peak1, peak2)
	gc()
	
	# check bmat	
	bmat1 = obj1@bmat;
	bmat2 = obj2@bmat;
	if((length(bmat1) == 0) != (length(bmat2) == 0)){
		stop("bmat has different dimentions in obj1 and obj2!")
	}else{
		bmat = base::rbind(bmat1, bmat2);
	}
	rm(bmat1, bmat2)
	gc()
	
	# check gmat	
	gmat1 = obj1@gmat;
	gmat2 = obj2@gmat;
	if((length(gmat1) == 0) != (length(gmat2) == 0)){
		stop("gmat has different dimentions in obj1 and obj2!")
	}else{
		gmat = base::rbind(gmat1, gmat2);
	}
	rm(gmat1, gmat2)
	gc()

	# check pmat	
	pmat1 = obj1@pmat;
	pmat2 = obj2@pmat;
	if((length(pmat1) == 0) != (length(pmat2) == 0)){
		stop("pmat has different dimentions in obj1 and obj2!")
	}else{
		pmat = base::rbind(pmat1, pmat2);
	}
	rm(pmat1, pmat2)
	gc()


	# check gmat	
	dmat1 = obj1@smat@dmat;
	dmat2 = obj2@smat@dmat;
	
	if((length(dmat1) == 0) != (length(dmat2) == 0)){
		stop("dmat has different dimentions in obj1 and obj2!")
	}else{
		dmat = base::rbind(dmat1, dmat2);
	}
	rm(dmat1, dmat2)
	gc()

	res = newSnap();
	res@feature = feature;
	res@barcode = c(obj1@barcode, obj2@barcode);
	res@file = c(obj1@file, obj2@file);
	res@sample = c(obj1@sample, obj2@sample);
	res@metaData = metaData;
	res@bmat = bmat;
	res@pmat = pmat;
	res@peak = peak;
	res@gmat = gmat;
	res@smat@dmat = dmat;
	res@smat@sdev = obj1@smat@sdev;
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
#' @examples
#' data(demo.sp);
#' filterCells(
#'	obj=demo.sp, 
#'	subset.names=c("UMI"), 
#'	low.thresholds=c(10),
#'	high.thresholds=c(Inf)
#'	);
#'
#' @return Returns a snap object containing only the relevant subset of cells
#' 
#' @importFrom methods is
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
		if(!is(obj, "snap")){
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

extractReadsFromOneCell <- function(
	barcode, 
	file
){
	if(missing(file)){
		stop("file is missing");
	}else{
		if(!isSnapFile(file)){
			stop(paste0(file, " is not a snap file"));			
		}
	}

	# read the reference barcode list
	barcode.list = as.character(tryCatch(barcode.list <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @extractReadsFromOneCell: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
	if(missing(barcode)){
		stop("barcode is missing");
	}else{
		if(class(barcode) != "character"){
			stop(paste0("barocde ", barcode, " is not character object"));
		}else{
			if(!(barcode %in% barcode.list)){
				stop(paste0("barocde ", barcode, " does not exist in snape file", file));
			}
		}		
	}
	
	pos.list = as.numeric(tryCatch(pos.list <- h5read(file, "FM/barcodePos"), error = function(e) {print(paste("Warning @readMetaData: 'FM/barcodePos' not found in ",file)); return(c())}));
	len.list = as.numeric(tryCatch(len.list <- h5read(file, "FM/barcodeLen"), error = function(e) {print(paste("Warning @readMetaData: 'FM/barcodeLen' not found in ",file)); return(c())}));
	
	if((x=length(pos.list)) != (y=length(len.list))){
		stop("FM/barcodeLen has different length with FM/barcodePos")
	}
	
	pos = pos.list[match(barcode, barcode.list)];
	len = len.list[match(barcode, barcode.list)];
	idx.arr = seq(pos, pos + len - 1);

	chroms = h5read(file, "FM/fragChrom", index=list(idx.arr))
	starts = h5read(file, "FM/fragStart", index=list(idx.arr))
	lens = h5read(file, "FM/fragLen", index=list(idx.arr))
	frags.gr = GRanges(chroms, 
		IRanges(starts, starts + lens - 1), 
		barcode=rep(barcode, length(chroms)), 
		file=rep(file, length(chroms))
		);
	return(frags.gr)
}

#' Extract Reads By Barcodes
#'
#' This function takes a barcode list and snap file as input and quickly extract reads belonging to the given barcodes. 
#' 
#' @param barcodes A vector contains the selected barcodes.
#' @param files A vector contains the snap file that barcodes belong to.
#' @param do.par A logic variable indicates weather to run this in parallel with multiple processors.
#' @param num.cores Number of processors to use.
#' 
#' @return Returns A GenomicRanges object that contains the reads
#' 
#' @importFrom methods is
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm
#' @importFrom methods slot
#' @export
extractReads <- function(
	barcodes,
	files,
	do.par=TRUE,
	num.cores=1
){
	if(missing(barcodes)){
		stop("barcodes is missing")
	}else{
		if(class(barcodes) != "character"){
			stop("barcodes must be character");
		}
	}

	if(missing(files)){
		stop("files is missing")
	}else{
		if(class(files) != "character"){
			stop("files must be character");
		}
	}
	ncell = length(barcodes);
	nfile = length(files);
	fileList = as.list(unique(files));
	if(ncell != nfile){
		stop("barcodes have different length with barcodes");
	}
	
	# check if snap files exist
	if(any(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){file.exists(x)})) == FALSE)
		print("error: these files does not exist")
		print(fileList[idx])
		stop()
	}
	
	# check if files are all snap files
	if(any(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){isSnapFile(x)})) == FALSE)
		print("error: these files are not snap file")
		print(fileList[idx])
		stop()
	}
	
	# check if FM session exist
	if(any(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)){
		idx = which(do.call(c, lapply(fileList, function(x){ "FM" %in% h5ls(x, recursive=1)$name  })) == FALSE)
		print("error: the following nsap files do not contain FM session")
		print(fileList[idx])
		stop()
	}

	if(do.par){
	    # input checking for parallel options
		if(num.cores > 1){
	        if (num.cores == 1) {
	          num.cores = 1
	        } else if (num.cores > detectCores()) {
	          num.cores <- detectCores() - 1
	          warning(paste0("num.cores set greater than number of available cores(", parallel::detectCores(), "). Setting num.cores to ", num.cores, "."))
	        }
	      } else if (num.cores != 1) {
	        num.cores <- 1
		}
	}
	read.ls = mclapply(as.list(seq(barcodes)), function(i){
		extractReadsFromOneCell(
			barcode = barcodes[i], 
			file = files[i]
		)		
	}, mc.cores=num.cores);
	read.gr = Reduce(c, read.ls);
	return(read.gr);
}

#' Feature filtration
#'
#' This function takes a snap obj as input and perform feature selection in the following manner:
#' 1) calculate coverage of each genomic bin/feature;
#' 2) log scale the coverage by log10(count + 1);
#' 3) the log-nromal distribution is then converted into zscore; 
#' 4) bins with zscore beyond [low.threshold, high.threshold] are filtered;
#' 
#' @param obj A snap obj
#' @param low.threshold Low cutoffs for the parameters (default is -2)
#' @param high.threshold High cutoffs for the parameters (default is 2)
#' @param mat Matrix to filter c("bmat", "pmat")
#'
#' @examples
#' data(demo.sp);
#' filterBins(
#'	obj=demo.sp, 
#'	low.threshold=-2,
#'	high.threshold=2
#'	);
#' 
#' @return Returns a snap obj
#' @importFrom stats sd
#' @importFrom methods slot is
#' @export
filterBins <- function(obj, low.threshold, high.threshold, mat) {
  UseMethod("filterBins", obj);
}

#' @export
filterBins.default <- function(obj, low.threshold=-2, high.threshold=2, mat=c("bmat", "pmat")){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("'obj' is not a snap obj")
		};		
	}
	
	mat = match.arg(mat);
	data.use = methods::slot(obj, mat);

	if((x=nrow(data.use)) == 0L){		
		stop("count matrix is empty")
	}

	idy = which(Matrix::colSums(data.use) > 0);
	cov = log(Matrix::colSums(data.use)[idy] + 1, 10);
	zcov = (cov - mean(cov)) / stats::sd(cov);	
	idy2 = which(zcov >= low.threshold & zcov <= high.threshold);
	idy = idy[idy2];
	methods::slot(obj, mat) = data.use[,idy,drop=FALSE];
	if(mat=="bmat"){
		obj@feature = obj@feature[idy];
	}else if(mat=="pmat"){
		obj@peak = obj@peak[idy];		
	}
	return(obj)
}

#' Check a snap-format file
#' 
#' This function takes a file name as input and check if the file is a snap-formated file 
#' @param file A file name.
#' @examples
#' file.name = system.file("extdata", "demo.snap", package = "SnapATAC");
#' isSnapFile(file.name);
#' 
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

readBins <- function(file, bin.size=5000){	
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
	if(!(bin.size %in% showBinSizes(file))){stop(paste("Error @addBmat: bin.size ", bin.size, " does not exist in ", file, "\n", sep=""))};
	options(scipen=999);
	binChrom = tryCatch(binChrom <- h5read(file, paste("AM", bin.size, "binChrom", sep="/")), error = function(e) {stop(paste("Warning @readaddBmatSnap: 'AM/bin.size/binChrom' not found in ",file))})
	binStart = tryCatch(binStart <- h5read(file, paste("AM", bin.size, "binStart", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})
	if(bin.size == 0){
		binEnd = tryCatch(binEnd <- h5read(file, paste("AM", bin.size, "binEnd", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/binStart' not found in ",file))})		
	}else{
		binEnd = binStart + as.numeric(bin.size) -1;
	}
	if((length(binChrom) == 0) || (length(binStart) == 0)){stop("Error @addBmat: bin is empty! Does not support empty snap file")}
	if(length(binChrom) != length(binStart)){
		stop(paste("Error @addBmat: ", "binChrom and binStart has different length!", sep=""))
	}else{
		nBin = length(binChrom);
	}
	bins = GRanges(binChrom, IRanges(as.numeric(binStart),binEnd), name=paste(paste(binChrom, binStart, sep=":"), binEnd, sep="-"));				
	rm(binChrom, binStart);
	return(bins)
}

readPeaks <- function(file){	
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	
	if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
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
	return(bins)
}

#' @importFrom methods is
#' @import Matrix
addBmatToSnapSingle <- function(obj, file, bin.size=5000){	
	# close the previously opened H5 file
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
	}

	if(!file.exists(file)){stop(paste("Error @addBmat: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @addBmat: ", file, " is not a snap-format file!", sep=""))};
	if(!(bin.size %in% showBinSizes(file))){stop(paste("Error @addBmat: bin.size ", bin.size, " does not exist in ", file, "\n", sep=""))};
	obj@bmat = Matrix(0,0,0);
	############################################################################
	barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @addBmat: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));	

	bin.sizeList = showBinSizes(file);
	if(length(bin.sizeList) == 0){stop("Error @addBmat: bin.sizeList is empty! Does not support reading empty snap file")}
	if(!(bin.size %in% bin.sizeList)){stop(paste("Error @addBmat: ", bin.size, " does not exist in bin.sizeList, valid bin.size includes ", toString(bin.sizeList), "\n", sep=""))}
	
	bins = readBins(file, bin.size);
	obj@feature = bins;
	idx = as.numeric(tryCatch(idx <- h5read(file, paste("AM", bin.size, "idx", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idx' not found in ",file))}));
	idy = as.numeric(tryCatch(idy <- h5read(file, paste("AM", bin.size, "idy", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/idy' not found in ",file))}));
	count = as.numeric(tryCatch(count <- h5read(file, paste("AM", bin.size, "count", sep="/")), error = function(e) {stop(paste("Warning @addBmat: 'AM/bin.size/count' not found in ",file))}));	

	if(!all(sapply(list(length(idx),length(idy),length(count)), function(x) x == length(count)))){stop("Error: idx, idy and count has different length in the snap file")}	
	nBarcode = length(barcode);
	nBin = length(obj@feature);
	M = sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nBin));
	rownames(M) = barcode;
	obj@bmat = M[match(obj@barcode, rownames(M)),]
	rm(idx, idy, count, M);
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	gc();
	return(obj);
}

#' @importFrom methods is
addGmatToSnapSingle <- function(obj, file){
        # close the previously opened H5 file
		if(exists('h5closeAll', where='package:rhdf5', mode='function')){
			rhdf5::h5closeAll();		
		}else{
			rhdf5::H5close();
		}
        # check the input
        if(!is(obj, "snap")){stop(paste("Error @addGmat: ", file, " does not exist!", sep=""))}
        if(!file.exists(file)){stop(paste("Error @addGmat: ", file, " does not exist!", sep=""))};
        if(!isSnapFile(file)){stop(paste("Error @addGmat: ", file, " is not a snap-format file!", sep=""))};

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
		nBarcode = length(barcode);
		nGene = length(geneName);
		M = sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nGene));
		rownames(M) = barcode;
		obj@gmat = 	M[match(obj@barcode, rownames(M)),]
		colnames(obj@gmat) = geneName;
		rm(idx, idy, count, M);

		if(exists('h5closeAll', where='package:rhdf5', mode='function')){
			rhdf5::h5closeAll();		
		}else{
			rhdf5::H5close();
		}
		gc();
		return(obj);
}

#' @importFrom methods is
addPmatToSnapSingle <- function(obj, file){
        # close the previously opened H5 file
		if(exists('h5closeAll', where='package:rhdf5', mode='function')){
			rhdf5::h5closeAll();		
		}else{
			rhdf5::H5close();
		}
        # check the input
        if(!is(obj, "snap")){stop(paste("Error @addPmat: ", file, " does not exist!", sep=""))}
        if(!file.exists(file)){stop(paste("Error @addPmat: ", file, " does not exist!", sep=""))};
        if(!isSnapFile(file)){stop(paste("Error @addPmat: ", file, " is not a snap-format file!", sep=""))};

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
		nBarcode = length(barcode);
		nPeak = length(obj@peak);
		M = sparseMatrix(i=idx, j =idy, x=count, dims=c(nBarcode, nPeak));
		rownames(M) = barcode;
		obj@pmat = 	M[match(obj@barcode, rownames(M)),]
		rm(idx, idy, count, M);
		if(exists('h5closeAll', where='package:rhdf5', mode='function')){
			rhdf5::h5closeAll();		
		}else{
			rhdf5::H5close();
		}
		
		gc();
		return(obj);
}

#' @importFrom methods is
createSnapSingle <- function(file, sample, metaData=TRUE, description=NULL){	
	# close the previously opened H5 file
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	# check the input
	if(!file.exists(file)){stop(paste(file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste(file, " is not a snap-format file!", sep=""))};
	if(!(is.logical(metaData))){stop(paste("metaData is not a logical variable!", sep=""))};
	
	if(!(is.null(description))){
		if(!is(description, "character")){
			stop("description must be character object")
		}
	}else{
		description=character()
	}
	
	# create an empty snap object
	res = newSnap();
	############################################################################
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
	res@sample = rep(sample, length(barcode));
	res@file = rep(normalizePath(file), length(res@barcode));
	if(exists('h5closeAll', where='package:rhdf5', mode='function')){
		rhdf5::h5closeAll();		
	}else{
		rhdf5::H5close();
	}
	gc();
	return(res);	
}

