#' Identifiy Differential Accessible Regions (DARs)
#'
#' This function takes a snap object and finds differentially 
#' accessible regions (DARs) that define clusters. 
#'
#' @param obj A snap object.
#' @param input.mat Matrix to use for finding differential features c("bmat", "pmat", "gmat").
#' @param cluster.pos Cluster to identify DAR markers.
#' @param cluster.neg Cluster used as negative control compare with cluster.pos [NULL].
#' If cluster.neg is NULL, runFindDARs will automatically identifies background cells 
#' by finding those that are closest to cluster.pos cells as a local background.
#' @param bcv Biological coefficient of variation. Typical values for the common BCV 
#' (square-rootdispersion) for datasets arising from well-controlled experiments are 
#' 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 
#' for technical replicates.
#' @param cluster.neg.method Method to find negative control cells if cluster.neg==NULL. ["knn", "random"]
#' @param test.method Test method for differential analysis c("exactTest", "LRT", "QLF").
#' @param seed.use Random seeds.
#' 
#' @examples
#' data(demo.sp);
#' idy = findDAR(
#'	obj=demo.sp, 
#'  input.mat="pmat", 
#'  cluster.pos=1, 
#'  bcv=0.1, 
#'  test.method="exactTest", 
#'  seed.use=10
#'  );
#' 
#' @importFrom methods slot
#' @importFrom edgeR DGEList glmLRT glmFit glmQLFit exactTest
#' @importFrom stats model.matrix
#' @export

findDAR <- function(obj, input.mat, cluster.pos, cluster.neg, cluster.neg.method, bcv, fdr, pvalue, test.method, seed.use) {
  UseMethod("findDAR", obj);
}

#' @export
findDAR.default <- function(
	obj,
	input.mat=c("pmat", "bmat", "gmat"),
	cluster.pos,
	cluster.neg=NULL,
	cluster.neg.method=c("knn", "random"),
	bcv=0.1,
	test.method=c("exactTest", "LRT", "QLF"),
	seed.use=10
	){
		message("Epoch: checking inputs ...")		
		if(missing(obj)){
			stop("obj is missing");
		}else{
			if(!is(obj, "snap")){
				stop("obj is not a snap obj")
			}
			if((x=length(levels(obj@cluster))) == 0L){
				stop("obj does not have cluster, 'runCluster' first");
			}
		}
		ncell = nrow(obj);

		cluster.neg.method = match.arg(cluster.neg.method);
		
		input.mat = match.arg(input.mat);
		mat.use = methods::slot(obj, input.mat);
		if((x=nrow(mat.use)) == 0L){
			stop("input.mat is empty, add matrix to snap object first")
		}
		
		if((x=nrow(mat.use)) != ncell){
			stop("input.mat has different length with cell number, re-add this matrix to snap object");
		}
		
		if(missing(cluster.pos)){
			stop("cluster.pos is missing")
		}else{
			if(any(!cluster.pos %in% levels(obj@cluster))){
				stop("cluster.pos does not exist in cluster")
			}
			idx.pos = which(obj@cluster %in% cluster.pos);
		}
		
		if(is.null(cluster.neg)){
			if(isKgraphEmpty(obj@graph)){
				stop("knn graph is empty, run runKNN first or reset cluster.neg")
			}
			idx.neg = findNegCells(obj, idx.pos, method=cluster.neg.method);
		}else{
			if(!(cluster.neg %in% levels(obj@cluster))){
				stop("cluster.neg does not exist in cluster")
			}
			idx.neg = which(obj@cluster == cluster.neg);
		}
		
		ncell.pos = length(idx.pos);
		ncell.neg = length(idx.neg);
		
		message("Epoch: identifying DARs for positive cluster ...")
		test.method = match.arg(test.method);
		
		mat.use.pos = mat.use[idx.pos,,drop=FALSE];
		mat.use.neg = mat.use[idx.neg,,drop=FALSE];
		
		# perform test
		data.use = data.frame(Matrix::colSums(mat.use.neg), Matrix::colSums(mat.use.pos))
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
		return(tb.pos);
	}
