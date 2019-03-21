#' @importFrom methods setClass setMethod
NULL
#' An S4 class jaccard to represent a knn graph object.
#'
#' Class kgraph defines a kgraph object.
#'
#' @slot mat a parase matrix object that contains the knn adjacent matrix
#' @slot k k used for constructing knn graph 
#' @slot file file name that saves that knn graph
#' @slot snn knn is converted to snn graph if True
#' @slot snn.prune snn prunning. edges with weight less than snn.prune will be removed
#' @name kgraph-class
#' @rdname kgraph-class
#' @importFrom methods setClassUnion
methods::setClassUnion("MatrixOrmatrix", c("Matrix", "matrix"))
kgraph <- setClass(
  Class = "kgraph",
  slots = list(
    mat = "MatrixOrmatrix",
	file = "character",
    k = "numeric",
	snn = "logical",
	snn.prune = "numeric"
  )
)

setMethod(
  f = 'show',
  signature = 'kgraph',
  definition = function(object) {
    cat(
      'number of cells:', nrow(object@mat), '\n',
      'K: ', object@k, '\n',
      'graph file: ', object@file, '\n',
      'snn: ', object@snn, '\n',
      'snn.prune: ', object@snn.prune, '\n'
    )
  }
)

#' subsetting for kgraph objects
#'
#' This function takes a kgraph object and returns the subset of kgraph object.
#' @param x A kgraph object
#' @param i selected rows
#' @param j selected columns
#' @param drop drop unused levels
#' @export
setMethod(
	f="[", 
    signature = 'kgraph',
	function(x,i,j, drop="missing"){
		.mat = x@mat;
		# a single row or column
       if(!missing(i)){
		   if(max(i) > nrow(.mat)){
			   stop("idx exceeds number of cells");
		   }
		   if(nrow(.mat) > 0){.mat <- .mat[i,i,drop=FALSE]}
	   }
	   if(!missing(j)){
		   stop("kgraph does not support subsetting for columns");
	   }
	   x@mat = .mat;
	   return(x);
})

#' @importFrom methods new
newKgraph <- function (mat=NULL, file=NULL, k=NULL, snn=NULL, snn.prune=NULL) {
	if(is.null(mat)){
		mat = Matrix::Matrix(0,0,0, sparse=TRUE);
	}
	if(is.null(file)){
		file = character();
	}

	if(is.null(k)){
		k = numeric();
	}

	if(is.null(snn)){
		snn = FALSE;
	}

	if(is.null(snn.prune)){
		snn.prune = numeric();
	}
	
	res = new("kgraph", 
			  mat=mat,
			  file=file,
			  k=k,
			  snn=snn,
			  snn.prune=snn.prune
			  );
	return(res)	
}


#' @importFrom methods new
isKgraphEmpty <- function(obj){
	if(is.null(obj)){
		stop("obj is empty")
	}else{
		if(class(obj) != "kgraph"){
			stop("obj is not a kgraph object");
		}else{
			if((x = nrow(obj@mat)) > 0L){
				return(FALSE)
			}
			if((x=length(obj@file) > 0L)){
				if(file.exists(obj@file)){
					return(FALSE)
				}
			}
		}
	}
	return(TRUE);
}

getGraph <- function(obj){
	if(is.null(obj)){
		stop("obj is empty")
	}else{
		if(class(obj) != "kgraph"){
			stop("obj is not a kgraph object");
		}
		if(isKgraphEmpty(obj)){
			stop("obj is empty");			
		}
	}
	
	if((x=nrow(obj@mat) != 0L)){
		return(obj@mat);
	}
	
	if(file.exists(obj@file)){
		edgeList = read.table(obj@file, header=FALSE);
		if((x=ncol(edgeList)) != 3L){
			stop(paste(obj@file, " does not have 3 columns"));
		}
		num.node = max(edgeList[,c(1,2)]);
		M1 = sparseMatrix(i=edgeList[,1], j=edgeList[,2], x=edgeList[,3], dims=c(num.node,num.node));
		M2 = sparseMatrix(i=edgeList[,2], j=edgeList[,1], x=edgeList[,3], dims=c(num.node,num.node));
		M = M1 + M2;
		rm(M1, M2);
		gc();
		return(M);
	}
}

# Write edge list to a file
# @param edges A data.frame obj contains edges and edge weight
# @param file Output file name
# @param ... Arguments passed to write.table
#' @importFrom utils write.table
writeEdgeListToFile <- function(edges, file, ...
){
	if(missing(edges) || missing(file)){
		stop("missing edges or file inputs")
	}
	
	if(!file.create(file)){
		stop("fail to create file")
	}			
	
    write.table(edges, file = file, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE)
}

readEdgeListFromFile <- function(file)
{
	edgeList = read.table(file, header=FALSE);
	return(edgeList);
}

