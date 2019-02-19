#' Frequency–inverse document frequency (IF-IDF)
#'
#' This function takes a Snap object as input and perform IF-IDF against the binary matrix.
#' 
#' Latent Semantic Indexing (LSI), a technique commonly used in the natural language processing 
#' that analyze the relationships between a set of documents, has been applied to scale the 
#' binary matrix and infer the association between cells. LSI first created a TF-ITD 
#' (frequency–inverse document frequency) matrix which scales the binary matrix by weighting 
#' each site accessible in an individual cell by the total number of sites accessible in that cell. 
#' 
#' @param object A snap object
#'
#' @export
#'
normTFIDF <- function(object, ...) {
  UseMethod("normTFIDF");
}

normTFIDF.default <- function(object){	
	
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("Not a snap object")
	}
		
	# 2. check of bmat exists;
	if(nrow(object@bmat) == 0){
		stop("@bmat is empty")		
	}	

	# 3. check of nmat exists;
	if(nrow(object@nmat) != 0){
		message("Warning: @nmat already exists")		
	}	
	
	x = Matrix::t(object@bmat);
	nfreqs = t(t(x) / Matrix::colSums(x));
	idf = as(log(1 + ncol(x) / Matrix::rowSums(x)), "sparseVector");
	tf_idf_counts = as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs;
	object@nmat = t(tf_idf_counts);	
}

# library(TailRank)
# tf_idf_counts
# x = tf_idf_counts
# x_mu = Matrix::rowMeans(x)
# x_var = c(apply(x[1:100000,], 1, sd), apply(x[100001:200000,], 1, sd), apply(x[200001:300000,], 1, sd), apply(x[300001:nrow(x),], 1, sd))
# x_cv = x_var/x_mu * 100
# linearModelVar <- lm(y ~ x, data.frame(x=log(x_mu), y=log(x_cv)))
# idy = which(linearModelVar$residuals > 0.5)


#RowVar <- function(x, ...) {
#  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
#}
