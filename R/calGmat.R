#' Calculate Cell x Gene Table from a snap Object
#'
#' This function takes a snap object and a GenomicRanges object that indicates the 
#' features or genes locations of interest and returns a cell x gene count table. 
#'  
#' @param object A Snap object.
#' @param gene A GenomicRange object that contains the genomic intervals of interest.
#' @param norm A logical class that indicates whether to normalize the count matrix. If True, the count table will be normalized to FPKM (default TRUE).
#' @param ncores A numeric class that indicates how many processers to use for computation (default 1). The more used, the faster it runs.
#'
#' @return Returns a Snap object with the cell x gene matrix stored in object@gmat
#'
#' @export
calCellGeneTable <- function(object, gene, ...) {
  UseMethod("calCellGeneTable", object);
}

#' @export
calCellGeneTable.default <- function(
	object, 
	gene, 
	ncores=1, 
	mat=c("bmat", "pmat")
){
	if(!is.snap(object)){
		stop("'object' is not a snap object")
	}

	if(class(gene) != "GRanges"){
		stop("'gene' is not a GenomicRanges object")
	}
		
	mat = match.arg(mat);
	if(mat == "bmat"){
		if(nrow(object@bmat) == 0){
			stop("'bmat' is empty")
		}		
		if(max(object@bmat) > 1){
			stop("'bmat' is not a binary matrix, run 'makeBinary' first")
		}
		feature = object@feature;
		mat = object@bmat;
	}else{
		if(nrow(object@pmat) == 0){
			stop("'pmat' is empty")
		}		
		feature = object@peak;
		mat = object@pmat;
	}
	
	if(any(duplicated(gene$gene_name))){
		stop("'gene' contains duplicate gene names, remove duplicate first");
	}
	
	# find overlap between feature and genes
	ov = data.frame(findOverlaps(feature, gene));
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
		}, mc.cores=ncores);
		# combine vectors and create a count matrix;
		count.mt = t(as.matrix(do.call(rbind, count.ls)));
		count_table = Matrix(0, nrow=nrow(object), ncol=length(gene), sparse=TRUE);
		count_table[,as.numeric(names(ov.ls))] = count.mt;
		count_table = count_table;		
	}else{
		count_table = Matrix(0, nrow(object), ncol=length(gene), sparse=TRUE);
	}
	colnames(count_table) = gene$name;
	object@gmat = count_table;
	colnames(object@gmat) = gene$gene_name;
	return(object);
}

