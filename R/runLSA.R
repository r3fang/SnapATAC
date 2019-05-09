#' Latent Semantic Analysis
#'
#' This function takes a Snap obj as input with bmat/pmat/gmat slot and run 
#' Latent Semantic Analysis (LSA).  
#' 
#' Below instruction is modified from 10X cell-ranger website
#' The a cell-by-bin (bmat) or cell-by-peak (pmat) matrix is first normalized via the inverse-document 
#' frequency (idf) transform where each peak/bin count is scaled by the log of the ratio of the number 
#' of barcodes in the matrix and the number of barcodes where the peak has a non-zero count. This provides 
#' greater weight to counts in peaks that occur in fewer barcodes. Singular value decomposition (SVD) is 
#' performed on this normalized matrix using IRLBA without scaling or centering, to produce the transformed 
#' matrix in lower dimensional space, as well as the components and the singular values signifying the 
#' importance of each component. 
#'
#' LSA has four major steps:
#' 1) term frequency - TF = t(t(X) / Matrix::colSums(X)); When logTF is TRUE, TF is also log scaled.
#' 2) inverse document frequency - IDF = log(1 + ncol(X) / rowSums(X))
#' 3) TF-IDF - TF * IDF
#' 4) SVD - Run singular value decomposition
#' 
#' @param obj A snap obj
#' @param input.mat Input matrix to be used for LSA c("bmat", "pmat").
#' @param pc.num An integer number of dimetions to return [50].
#' @param logTF A logical variable indicates wehther to log-scale term frequency [TRUE].
#' @param scale.factor A numeric variable used to scale the logTF [100000].
#' @param seed.use A numeric class that indicates random seeding number [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp);
#' demo.sp = runLSA(obj=demo.sp, input.mat="bmat", pc.num=50, logTF=TRUE, min.cell=0);
#' 
#' @importFrom irlba irlba
#' @import Matrix
#' @export

runLSA <- function(
	obj, 
	input.mat=c("bmat", "pmat"), 
	pc.num=50, 
	logTF=TRUE, 
	scale.factor=100000,
	min.cell=10,
	seed.use=10
){
	set.seed(seed.use);
	message("Epoch: checking the inputs ...");
	
	# 1. check if obj is a snap;
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
	}
	
	if(!is(pc.num, "numeric")){
		stop("pc.num is not a numeric")
	}
	
	if(!is(logTF, "logical")){
		stop("logTF is not a numeric")
	}
	
	# 2. check if input matrix exists
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use = t(obj@bmat);
	}else if(input.mat == "pmat"){
		data.use = t(obj@pmat);
	}else if(input.mat == "gmat"){
		data.use = t(obj@gmat);
	}else{
		stop("input.mat is not supported!")
	}
	
	if((x = max(data.use@x)) > 1L){
		stop("input matrix is not a binary matrix, run 'makeBinary' first");
	}
	
	ncell = ncol(data.use);
	if(!is(min.cell, "numeric")){
		stop("min.cell is not a numeric")
	}else{
		if(min.cell > ncell){
			stop("min.cell exceeds number of cells")
		}
	}
	
	message("Epoch: filtering low-coverage matrix ...");
	idx = which(Matrix::rowSums(data.use) > min.cell);
	if((x = length(idx)) == 0L){
		stop("No features has coverage than min.cell")
	}else{
		data.use = data.use[idx,];
	}
	
	cell.cov = Matrix::colSums(data.use);
	if(any(cell.cov == 0)){
		stop("empty cells detected after removing low-coverage features, remove empty cells first.")
	}
	
	message("Epoch: running term frequency ...");		
	tf = t(t(data.use) / Matrix::colSums(data.use));
	
	if(logTF){
		message("Epoch: running log term frequency ...");
        tf@x = log1p(tf@x * scale.factor);
	}
	
	message("Epoch: running inverse document frequency ...");
    idf = log(1 + ncol(data.use) / Matrix::rowSums(data.use))

	message("Epoch: running TF-IDF ...");
    tf = t(tf);
    tf@x <- tf@x * rep.int(idf, diff(tf@p));
    tf = t(tf);
	
	message("Epoch: running SVD ...");
	set.seed(seed.use);
	pca.results = irlba(t(tf), nv=pc.num);
	obj@smat@dmat = pca.results$u;
	obj@smat@sdev = pca.results$d;
	obj@smat@iter = pca.results$iter;
	obj@smat@imat = input.mat;
	
	if(logTF){
		obj@smat@method = "LSA.logTF";				
	}else{
		obj@smat@method = "LSA";		
	}
	message("Epoch: Done");		
	return(obj);
}
