#' Identifiy Statistically Singficant Differential Accessible Regions (DARs)
#'
#' This function takes a snap object and a cluster vector and replicate indicator as input
#' and identify the differential accessibile regions (DARs) features between 2 groups of cells.
#' 
#' @param object A Snap object
#' @param cluster A vector that indicates the cluster label for each cell
#' @param rep A vector that indicates the replicate label for each cell
#' @param treat A character or numeric object indicates which cluster or clusters are used as treatment group
#' @param ctrl A character or numeric object indicates which cluster or clusters are used as control group
#' @param method A character object indicates what test to use (LRT, EXACT, QLF)
#' @param down.sample A logic object indicates whether or not to perform down sampling of bigger cluster
#' @param dispersion A numberic object indicates the common dispersion for EXACT test
#' @param min.row.sum A numberic object that filters any feature with coverage less than min.row.sum
#' @param max.row.sum A numberic object that filters any feature with coverage more than min.row.sum
#' @param rand.seed A numeric object that sets the random seeds
#'
#' @return Returns a data.frame table that contains differential analysis info for each feature:
#' id - index for original region
#' cpmT - cpmC - averaged CPM for treatment and control
#' logFC - log2 fold change between treatment and control
#' logCPM - log2 average CPM between treatment and control
#' LR - Likilihood Ratio
#' PValue - Chi-square P-value
#' fdr - Adjusted P-value
#'
#' @export

findDAR <- function(object, input=c("cmat", "gmat"), cluster=NULL, rep=NULL, treat=NULL, ctrl=NULL, down.sample=TRUE, dispersion=0.2, min.row.sum = 10, max.row.sum=500, rand.seed=10, ...) {
  UseMethod("findDAR");
}

findDAR.default <- function(object, input=c("cmat", "gmat"), cluster=NULL, rep=NULL, treat=NULL, ctrl=NULL, down.sample=TRUE, dispersion=0.3, min.row.sum=10, max.row.sum=500, rand.seed=10, ...){
	# 1. check if object is a snap;
	if(class(object) != "snap"){
		stop("Not a snap object")
	}
	if(missing(cluster)){
		stop("cluster is missing")		
	}

	if(class(cluster) != "factor"){
		stop("cluster is not a factor")		
	}
	
	if(missing(rep)){
		stop("rep is missing")		
	}

	if(class(rep) != "factor"){
		stop("rep is not a factor")		
	}
	
	if((nrow(object) != length(cluster)) | nrow(object) != length(rep)){
		stop("cluster or rep incorrect length")				
	}

	if(any(treat %in% cluster == FALSE)){
		stop("treat is not contained in cluster")						
	}

	if(any(ctrl %in% cluster == FALSE)){
		stop("ctrl is not contained in cluster")						
	}
	
	input = match.arg(input);
	cluster = droplevels(cluster);
	rep = droplevels(rep);
	
	if(input == "cmat"){
		cmat = object@cmat;		
	}else if(input == "bmat"){
		cmat = object@bmat;		
	}else if(input == "gmat"){
		cmat = gmat;
	}
	
	nrep = length(levels(rep));
	if(nrep == 1){
		stop("Does not support 1 replciate!")
	}

	ncluster = length(levels(cluster));
	nfeature = ncol(cmat);
	
	if(nfeature <= 0 ){
		stop("0 feature detected in object")
	}

	# find the min cell number
	if(down.sample){
		ncell.treat <- min(length(which(cluster %in% treat)), length(which(cluster %in% ctrl)));
		ncell.ctrl <- min(length(which(cluster %in% treat)), length(which(cluster %in% ctrl)));
		set.seed(rand.seed);
		idx.treat <- sample(which(cluster %in% treat), ncell.treat);
		set.seed(rand.seed);
		idx.ctrl <- sample(which(cluster %in% ctrl), ncell.ctrl);
		# reconstruct all the input
	}else{
		ncell.treat <- length(which(cluster %in% treat));
		ncell.ctrl <- length(which(cluster %in% ctrl));
		set.seed(rand.seed);
		idx.treat <- sample(which(cluster %in% treat), ncell.treat);
		set.seed(rand.seed);
		idx.ctrl <- sample(which(cluster %in% ctrl), ncell.ctrl);
	}
		
	cmat = cmat[c(idx.treat, idx.ctrl),];
	cluster = cluster[c(idx.treat, idx.ctrl)];
	rep = rep[c(idx.treat, idx.ctrl)];
	
	# create the ensemble signal
	x = matrix(0, nfeature, 2*nrep)
	
	i = 1
	for(rep_i in levels(rep)){	
		x[,i] = Matrix::colSums(cmat[which(cluster %in% treat & rep == rep_i),])
		i = i + 1
	}

	for(rep_i in levels(rep)){		
		x[,i] = Matrix::colSums(cmat[which(cluster %in% ctrl & rep == rep_i),])
		i = i + 1
	}
	
	group <- factor(c(rep(1, nrep), rep(2, nrep)));
	design <- model.matrix(~group);
	keep = which(rowSums(x) >= min.row.sum & rowSums(x) <= max.row.sum)
	x = x[keep,,drop=FALSE]
	y <- DGEList(counts=x, group=group);	
	# pre-filteration of low-coverage and high coverage
	y <- calcNormFactors(y);
	
	y <- estimateDisp(y, design, trend.method="locfit");
	fit <- glmFit(y, design);				
	tb <- glmLRT(fit,coef=2)$table;
	
	#if(nrep == 1){
	#	if(method == "LRT"){
	#		fit <- glmFit(y, design, dispersion=0.2^2);				
	#		tb <- glmLRT(fit,coef=2)$table;
	#	}else if(method == "QLF"){		
	#		fit <- glmQLFit(y, design, dispersion=0.2^2)$table;				
	#	}else{
	#		tb <- exactTest(y, dispersion=0.2^2)$table;
	#	}	
	#}else{
	#	if(method == "LRT"){
	#		y <- estimateDisp(y, design, trend.method="locfit");
	#		fit <- glmFit(y, design);				
	#		tb <- glmLRT(fit,coef=2)$table;
	#	}else if(method == "QLF"){		
	#		y <- estimateDisp(y,design, trend.method="locfit");
	#		fit <- glmQLFit(y, design)$table;				
	#	}else{
	#		y <- estimateDisp(y,design, trend.method="locfit");
	#		tb <- exactTest(y)$table;
	#	}	
	#}
	tb$fdr <- p.adjust(tb$PValue, method="BH");
	return(cbind(id=keep, cpm(y), tb));
}
