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
  UseMethod("findDAR", object);
}

#' @export
findDAR.default <- function(
	object,
	mat=c("pmat", "bmat", "gmat"),
	barcodes.sel,
	bcv=0.1,
	rand.seed=10,
	pca_dims,
	fdr=5e-2,
	background_method=c("exact", "fast"),
	method=c("exactTest", "LRT", "QLF"),
	...
	){
		# 1. check if object is a snap;
		if(class(object) != "snap"){
			stop("object is not a snap object")
		}
		
		if(missing(barcodes.sel)){
			stop("barcodes.sel is missing")
		}else{
			if(any(!(barcodes.sel %in% object@barcode))){
				stop("some 'barcodes' do not exist in the object");
			}
		}
		
		mat = match.arg(mat);
		if(mat == "bmat"){
			cmat = object@bmat;
		}else if(mat == "pmat"){
			cmat = object@pmat;
		}else if(mat == "gmat"){
			cmat = object@gmat;
		}
		
		if(missing(pca_dims)){
			stop("pca_dims is missing")
		}
		
		method = match.arg(method);
		
		message("Identifying accessible regions using postive sample");
		# positive cells
		idx.pos = which(object@barcode %in% barcodes);
		
		# identify negative control cells
		background_method = match.arg(background_method);
		if(background_method == "exact"){
			d = as.matrix(dist(object@smat[,pca_dims]))
			idx.neg = order(colMeans(d[idx.pos,-idx.pos]))[1:min(ncol(d) - length(idx.pos), length(idx.pos))]
		}else if(background_method == "fast"){
			k = 50;
			dx = nn2(object@smat[,pca_dims], k = k+1)$nn.idx;
			dx = dx[,2:(k+1)];
			edges = matrix(unlist(sapply(1:nrow(dx),function(i) { rbind(rep(i,k),dx[i,])})),nrow=2);
			edges = as.matrix(t(edges));
			freq = table(edges[which((edges[,1] %in% idx.pos) & !(edges[,2] %in% idx.pos)),2]);
			idx.neg = as.numeric(names(freq)[order(freq, decreasing=TRUE)])[seq(min(length(idx.pos), length(freq)))];			
		}
		
		
		# calcualte coverage for posive and negative cells
		cmat.pos = cmat[idx.pos,];
		cmat.neg = cmat[idx.neg,];
		
		# perform test
		x = data.frame(colSums(cmat.neg), colSums(cmat.pos))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=x, group=group);
		
		if(method == "LRT"){
			fit <- glmFit(y, design, dispersion=bcv^2);
			tb.pos <- glmLRT(fit,coef=2)$table;
		}else if(method == "QLF"){
			tb.pos <- glmQLFit(y, design, dispersion=bcv^2)$table;
		}else{
			tb.pos <- exactTest(y, dispersion=bcv^2)$table;
		}
		
		message("Identifying accessible regions using negative control sample")
		# negative control by randomly select k cells
		set.seed(rand.seed);
		neg.idx.pos = sample(seq(nrow(object)), length(idx.pos))
		neg.idx.neg = sample(seq(nrow(object)), length(idx.pos))
				
		# calcualte coverage for posive and negative cells
		neg.cmat.pos = cmat[neg.idx.pos,];
		neg.cmat.neg = cmat[neg.idx.neg,];
		
		x = data.frame(colSums(neg.cmat.pos), colSums(neg.cmat.neg))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=x, group=group);
		
		if(method == "LRT"){
		fit <- glmFit(y, design, dispersion=bcv^2);
			tb.neg <- glmLRT(fit,coef=2)$table;
		}else if(method == "QLF"){
			tb.neg <- glmQLFit(y, design, dispersion=bcv^2)$table;
		}else{
			tb.neg <- exactTest(y, dispersion=bcv^2)$table;
		}
		
		message("calculating p-value and FDR table")
		fdr_table = data.frame()
		for(p_i in seq(0.001, 0.05, by=0.001)){
		fdr_table = rbind(fdr_table,
			data.frame(
			PValue=p_i,
			fdr=length(which(tb.neg$PValue < p_i & tb.neg$logFC > 0)) / length(which(tb.pos$PValue < p_i & tb.pos$logFC > 0))
			)
			)
		}
		if(length(which(fdr_table[,2] <= fdr)) > 0){
			fdr.idx = max(which(fdr_table[,2] <= fdr));
			return(which(tb.pos$PValue < fdr_table[fdr.idx,1] & tb.pos$logFC > 0))
		}
		return(c())
}