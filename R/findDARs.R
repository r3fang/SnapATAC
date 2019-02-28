#' Identifiy Statistically Singficant Differential Accessible Regions (DARs)
#' @export
findDAR <- function(object, ...) {
  UseMethod("findDAR", object);
}

#' @export
findDAR.default <- function(
	object,
	mat=c("pmat", "bmat", "gmat"),
	barcodes.sel,
	bcv=0.1,
	pca_dims=NULL,
	fdr=5e-2,
	background_method=c("knn", "random", "all", "other"),
	test_method=c("exactTest", "LRT", "QLF"),
	rand.seed=10,
	...
	){
		if(class(object) != "snap"){
			stop("object is not a snap object")
		}
		
		ncell = nrow(object);
		if(missing(barcodes.sel)){
			stop("barcodes.sel is missing")
		}else{
			if(any(!(barcodes.sel %in% object@barcode))){
				stop("some 'barcodes' do not exist in the object");
			}
		}
		
		mat = match.arg(mat);
		background_method = match.arg(background_method);
		test_method = match.arg(test_method);
		
		if(mat == "bmat"){
			cmat = object@bmat;
		}else if(mat == "pmat"){
			cmat = object@pmat;
		}else if(mat == "gmat"){
			cmat = object@gmat;
		}
		
		# check PCA dimentions
		if(is.null(pca_dims)){
			pca_dims = seq(ncol(object@smat));
		}else{
			if(any((pca_dims %in% seq(ncol(object@smat))) == FALSE)){
				stop("'pca_dims' exceeds the PCA space, reset 'pca_dims' and run it again")
			}		
		}
		
		message("Identifying accessible regions using postive sample");
		# positive cells
		idx.pos = which(object@barcode %in% barcodes.sel);
		
		# identify negative control cells
		# use exact method with cell number less than 5000
		if(background_method == "knn"){
			idx.neg = setdiff(seq(ncell), idx.pos)
			avg.profile = t(colMeans(object@smat[idx.pos,pca_dims]));
			nn.num = min(length(idx.pos), length(idx.neg));
			dx = nn2(object@smat[-idx.pos,pca_dims], avg.profile, nn.num)$nn.idx;
			idx.neg = idx.neg[dx[1,]];
		}else if(background_method == "random"){
			idx.neg = setdiff(seq(ncell), idx.pos);
			nn.num = min(length(idx.pos), length(idx.neg));
			idx.neg = sample(idx.neg, nn.num);	
		}else if(background_method == "all"){
			idx.neg = seq(ncell)
		}else if(background_method == "other"){
			idx.neg = setdiff(seq(ncell), idx.pos);			
		}
		
		# calcualte coverage for posive and negative cells
		cmat.pos = cmat[idx.pos,];
		cmat.neg = cmat[idx.neg,];
		
		# perform test
		x = data.frame(Matrix::colSums(cmat.neg), Matrix::colSums(cmat.pos))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=x, group=group);
		
		if(test_method == "LRT"){
			fit <- glmFit(y, design, dispersion=bcv^2);
			tb.pos <- glmLRT(fit,coef=2)$table;
		}else if(test_method == "QLF"){
			tb.pos <- glmQLFit(y, design, dispersion=bcv^2)$table;
		}else{
			tb.pos <- exactTest(y, dispersion=bcv^2)$table;
		}
		
		message("Identifying accessible regions using negative control sample")
		# negative control by randomly select k cells
		set.seed(rand.seed);
		neg.idx.pos = sample(seq(nrow(object)), length(idx.pos))
		background_method = "random";
		if(background_method == "knn"){
			neg.idx.neg = setdiff(seq(ncell), neg.idx.pos)
			avg.profile = t(colMeans(object@smat[neg.idx.pos,pca_dims]));
			nn.num = min(length(neg.idx.neg), length(neg.idx.pos));
			dx = nn2(object@smat[-neg.idx.pos,pca_dims], avg.profile, nn.num)$nn.idx;
			neg.idx.neg = neg.idx.neg[dx[1,]];
		}else if(background_method == "random"){
			neg.idx.neg = setdiff(seq(ncell), neg.idx.pos);
			nn.num = min(length(neg.idx.pos), length(neg.idx.neg));
			set.seed(rand.seed);
			neg.idx.neg = sample(neg.idx.neg, nn.num);	
		}else if(background_method == "all"){
			neg.idx.neg = seq(ncell)
		}else if(background_method == "other"){
			neg.idx.neg = setdiff(seq(ncell), neg.idx.pos);			
		}
						
		# calcualte coverage for posive and negative cells
		neg.cmat.pos = cmat[neg.idx.pos,];
		neg.cmat.neg = cmat[neg.idx.neg,];
		
		x = data.frame(Matrix::colSums(neg.cmat.pos), Matrix::colSums(neg.cmat.neg))
		group <- factor(c(1,2));
		design <- model.matrix(~group);
		y <- DGEList(counts=x, group=group);
		
		if(test_method == "LRT"){
			fit <- glmFit(y, design, dispersion=bcv^2);
			tb.neg <- glmLRT(fit,coef=2)$table;
		}else if(test_method == "QLF"){
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


