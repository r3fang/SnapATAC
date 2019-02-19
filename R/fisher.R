#' Identify Statistically Singficant Differential Accessible Regions Between 2 Clusters Using Fisher Exact Test
#'
#' This function takes a snap object and index (idx and idy) for control and treatment cells as input and identify differential accessibile regions (DARs) features between 2 groups or clusters of cells.
#' 
#' @param object A Snap object
#' @param input A character object that indicates what attribute of snap object to use for analysis ("cmat" by default, or "gmat", "bmat")
#' @param idx A vector that contains index for experiment group of cells
#' @param idy A vector that contains index for control group of cells
#' @param psudo.count A logical variable that indicates whether to add psudo count
#' @param min.cov A numeric variable between 0 to 1 that indicates min coverage (0.01) for feature to be tested. Any features are covered less than min.cov of all cells are not going to be tested
#' @param max.cov A numeric variable between 0 to 1 that indicates max coverage (0.6) for feature to be tested. Any features are covered more than max.cov of all cells are not going to be tested
#' @param min.fc A numeric variable that indicates min fold change (2) for feature to be tested. Any features with fold change less than min.fc will be not be tested
#' @param norm A logical variable that indicates whether to adjust for sequencing depth between two clusters
#' @param step A numeric variable used when norm is TRUE [1]
#'
#' @return Returns a data.frame table that contains differential analysis info for feature:
#' id - index for tested features (because not all features have been tested)
#' p1 - feature coverage for treatment group
#' p2 - feature coverage for control group (adjusted for sequencing depth if norm is TRUE)
#' p - average coverage among both control and treatment group (p1+p2)/2
#' log2FC - log2 fold change between treatment and control
#' pval - Fisher exact test P-value
#' @export

findDAR.fisher <- function(object, input=c("cmat", "gmat"), idx=NULL, idy=NULL, psudo.count=TRUE, min.cov=0.01, min.fc=2, max.cov = 0.8, norm=TRUE, step=1){
  UseMethod("findDAR.fisher", object);
}

#' @export
findDAR.fisher.default <- function(object, input=c("cmat", "gmat"), idx=NULL, idy=NULL, psudo.count=TRUE, min.cov=0.01, max.cov = 0.8, min.fc=2, norm=TRUE, step=1){
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("Not a snap object")
	}
    
	# 2. check if idx exceeds cell index; 
	if(any(!(idx %in% 1:nrow(object)))){
		stop("idx exceeds number of cells")
	}

	# 3. check if idy exceeds cell index; 
	if(any(!(idy %in% 1:nrow(object)))){
		stop("idy exceeds number of cells")
	}
	
	# 4. choose feature matrix to use
	input = match.arg(input);	
	if(input == "cmat"){
		cmat = object@cmat;		
	}else if(input == "bmat"){
		cmat = object@bmat;	
	}else if(input == "gmat"){
		cmat = object@gmat;
	}

	if(nrow(cmat) == 0){
		stop("Empty input matrix")
	}
	
	if(min.fc < 1){
		stop("min.fc has to be greater than 1");
	}
	
	nfeature = ncol(cmat);
	n1 = length(idx);
	n2 = length(idy);
	cmat1 = cmat[idx,]
	cmat2 = cmat[idy,]
	
	# count matrix binarization
	cmat1[Matrix::which(cmat1 > 0, arr.ind=TRUE)] = 1;
	cmat2[Matrix::which(cmat2 > 0, arr.ind=TRUE)] = 1;
	cmat = rbind(cmat1, cmat2)
		
	# remove features with very low coverage
	feature_sel = which(Matrix::colMeans(cmat) >= min.cov & Matrix::colMeans(cmat) <= max.cov)
	cmat1 = cmat1[,feature_sel]
	cmat2 = cmat2[,feature_sel]
	cmat = cmat[,feature_sel]
	
	count_matrix = data.frame(Matrix::colSums(cmat1)+1, n1+2, Matrix::colSums(cmat2)+1, n2+2)		

	# adjust for n2
	if(norm){
		message("Estimating delta")
		delta = 0
		p1 = count_matrix[,1]/count_matrix[,2]
		p2 = count_matrix[,3]/(count_matrix[,4] + delta)
		
		if(mean(log(p1/p2, 2)) < 0){
			while(mean(log(p1/p2, 2)) < 0){
				delta = delta + step
				p1 = count_matrix[,1]/count_matrix[,2]
				p2 = count_matrix[,3]/(count_matrix[,4] + delta)
			}
			
		}else{
			while(mean(log(p1/p2, 2)) >= 0){
				delta = delta - step
				p1 = count_matrix[,1]/count_matrix[,2]
				p2 = count_matrix[,3]/(count_matrix[,4] + delta)
			}
		}
		message(paste("Delta = ", delta, sep=""))
		#n2 = n2 + delta
		count_matrix[,4] = count_matrix[,4] + delta;
		#count_matrix = data.frame(colSums(cmat1)+1, n1+2, colSums(cmat2)+1, n2+2)		
	}

	count_matrix[count_matrix[,4] < count_matrix[,3],3] = count_matrix[count_matrix[,4] < count_matrix[,3],4];
	count_matrix$logFC = log(((count_matrix[,1])/(count_matrix[,2])) / ((count_matrix[,3])/(count_matrix[,4])), 2);
	idx = which(abs(count_matrix$logFC) >= log(min.fc, 2))
	feature_sel = feature_sel[idx]
	count_matrix = count_matrix[idx,]
	
	# calcualte pvalue
	pvalues = apply(count_matrix, 1, function(x) fisher.test(matrix(c(x[1], x[2] - x[1], x[3], x[4] - x[3]), 2, 2))$p.value)

	# adjust for fdr
	count_matrix$pvalues = pvalues
	count_matrix$p1 = count_matrix[,1]/count_matrix[,2];
	count_matrix$p2 = count_matrix[,3]/(count_matrix[,4]);
	count_matrix$p = (count_matrix$p1 + count_matrix$p2)/2;
	res = data.frame(id = feature_sel, p1 = count_matrix$p1, p2=count_matrix$p2, p=count_matrix$p, log2FC=count_matrix$logFC, pval=count_matrix$pvalue);
	return(res);
}


