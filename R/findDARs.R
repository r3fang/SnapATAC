#' Identifiy Differential Accessible Regions (DARs)
#'
#' This function takes a snap object and finds differentially 
#' accessible regions (DARs) that define clusters. 
#'
#' @param obj A snap object.
#' @param mat Matrix to use for finding differential features c("bmat", "pmat", "gmat").
#' @param cluster.pos Cluster to identify DAR markers.
#' @param cluster.neg Cluster used as negative control compare with cluster.pos [NULL].
#' If cluster.neg is NULL, runFindDARs will automatically identifies background cells 
#' by finding those that are closest to cluster.pos cells as a local background.
#' @param bcv Biological coefficient of variation. Typical values for the common BCV 
#' (square-rootdispersion) for datasets arising from well-controlled experiments are 
#' 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 
#' for technical replicates.
#' @param test.method Test method for differential analysis c("exactTest", "LRT", "QLF").
#' @param seed.use Random seeds.
#' @param fdr False discovery rate (FDR) [5e-2].
#' @param pvalue Pvalue [1e-2].
#' 
#' @importFrom methods slot
#' @importFrom edgeR DGEList glmLRT glmFit glmQLFit exactTest
#' @importFrom stats model.matrix
#' @export

findDAR <- function(obj, mat, cluster.pos, cluster.neg, bcv, fdr, pvalue, test.method, seed.use) {
  UseMethod("findDAR", obj);
}

#' @export
findDAR.default <- function(
	obj,
	mat=c("pmat", "bmat", "gmat"),
	cluster.pos,
	cluster.neg=NULL,
	bcv=0.1,
	fdr=5e-2,
	pvalue=5e-2,
	test.method=c("exactTest", "LRT", "QLF"),
	seed.use=10
	){
		message("Epoch: checking inputs ...")		
		if(missing(obj)){
			stop("obj is missing");
		}else{
			if(class(obj) != "snap"){
				stop("obj is not a snap obj")
			}
			if((x=length(levels(obj@cluster))) == 0L){
				stop("obj does not have cluster, 'runCluster' first");
			}
		}
		ncell = nrow(obj);
		
		mat = match.arg(mat);
		mat.use = methods::slot(obj, mat);
		if((x=nrow(mat.use)) == 0L){
			stop("mat is empty, add matrix to snap object first")
		}
		
		if((x=nrow(mat.use)) != ncell){
			stop("mat has different length with cell number, re-add this matrix to snap object");
		}
		
		if(missing(cluster.pos)){
			stop("cluster.pos is missing")
		}else{
			if(!(cluster.pos %in% levels(obj@cluster))){
				stop("cluster.pos does not exist in cluster")
			}
			idx.pos = which(obj@cluster == cluster.pos);
		}
		
		if(is.null(cluster.neg)){
			if(isKgraphEmpty(obj@kgraph)){
				stop("knn graph is empty, run runKNN first or reset cluster.neg")
			}
			idx.neg = findNegCells(obj, idx.pos, method="knn");
		}else{
			if(!(cluster.neg %in% levels(obj@cluster))){
				stop("cluster.neg does not exist in cluster")
			}
			idx.neg = which(obj@cluster == cluster.neg);
		}
		
		ncell.pos = length(idx.pos);
		ncell.neg = length(idx.neg);
		
		##############################################################################
		##############################################################################
		##############################################################################
		##############################################################################
		message("Epoch: identifying DARs for positive cluster ...")
		test.method = match.arg(test.method);
		
		mat.use.pos = mat.use[idx.pos,,drop=FALSE];
		mat.use.neg = mat.use[idx.neg,,drop=FALSE];
		
		# perform test
		data.use = data.frame(Matrix::colSums(mat.use.pos), Matrix::colSums(mat.use.neg))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=data.use, group=group);

		if(test.method == "LRT"){
			fit <- glmFit(y, design, dispersion=bcv^2);
			tb.pos <- glmLRT(fit,coef=2)$table;
		}else if(test.method == "QLF"){
			tb.pos <- glmQLFit(y, design, dispersion=bcv^2)$table;
		}else{
			tb.pos <- exactTest(y, dispersion=bcv^2)$table;
		}
		
		##############################################################################
		##############################################################################
		##############################################################################
		##############################################################################
		message("Epoch: identifying DARs for negative cluster ...")
		set.seed(seed.use);
		idx.pos = sample(seq(ncell), ncell.pos);
		idx.neg = findNegCells(obj, idx.pos, method="random");
		
		mat.use.pos = mat.use[idx.pos,,drop=FALSE];
		mat.use.neg = mat.use[idx.neg,,drop=FALSE];
		
		# perform test
		data.use = data.frame(Matrix::colSums(mat.use.pos), Matrix::colSums(mat.use.neg))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=data.use, group=group);

		if(test.method == "LRT"){
			fit <- glmFit(y, design, dispersion=bcv^2);
			tb.neg <- glmLRT(fit,coef=2)$table;
		}else if(test.method == "QLF"){
			tb.neg <- glmQLFit(y, design, dispersion=bcv^2)$table;
		}else{
			tb.neg <- exactTest(y, dispersion=bcv^2)$table;
		}
		
		message("Epoch: converting p-value to FDR ...")
		p1 = tb.pos$PValue;
		p1[which(tb.pos$logFC > 0)] = 1;

		p2 = tb.neg$PValue;
		p2[which(tb.neg$logFC > 0)] = 1;
		
		fdr_table = pvalue2fdr(p1, p2);
				
		if((x=length(which(fdr_table[,2] <= fdr))) > 0){
			fdr.idx = max(which(fdr_table[,2] <= fdr));
			pvalue = min(pvalue, fdr_table[fdr.idx,1]);
			return(which(p1 <= pvalue));
		}
		return(c())
}

