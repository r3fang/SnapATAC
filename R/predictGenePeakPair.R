#' Link Distal Regulatory Elements (peaks) to Putative Target Genes
#'
#' @param obj A snap object.
#' @param input.mat Input matrix c("pmat", "bmat").
#' @param gene.name Name of the target gene. It must be a gene in the cell-by-gene matrix [NULL].
#' @param gene.loci A genomic range object that contains the flanking window of the target gene [NULL].
#' @param do.par A logical variable indicates if to run this in parallel using multiple processors [FALSE].
#' @param num.cores A numeric class that indicates the number of cores to use for calculation [1].
#' @return Return A data.frame that contains the peaks and associating P-value.
#'
#' @importFrom stats glm
#' @import parallel
#' @import doSNOW
#' @importFrom plyr llply
#' 
#' @export
predictGenePeakPair <- function(
	obj, 
	input.mat=c("pmat", "bmat"),
	gene.name=NULL,
	gene.loci=NULL,
	do.par=FALSE,
	num.cores=1
){
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
	}
	
	
	# check gmat
	gene.mat = obj@gmat;
	if((x=nrow(gene.mat)) == 1L){
		stop("cell by gene matrix is empty!");
	}
	
	# check if gene.name exist in file
	if(!(gene.name %in% colnames(gene.mat))){
		stop("gene.name does not exist in the cell by gene matrix")
	}else{
		gene.val = gene.mat[,gene.name];
	}
	
	
	input.mat = match.arg(input.mat);	
	if(input.mat == "bmat"){
		data.use = obj@bmat;
		peak.use = obj@feature;
	}else if(input.mat == "pmat"){
		data.use = obj@pmat;
		peak.use = obj@peak;
	}else{
		stop("input.mat does not exist in obj")
	}
	
	if((x=max(data.use)) > 1L){
		stop("input matrix is not binarized, run 'makeBinary' first!")
	}

	if((x=length(peak.use)) == 1L){
		stop("peak is empty!")
	}

	if((x=nrow(data.use)) == 1L){
		stop("input matrix is empty!")
	}
	
	ncell = nrow(data.use);
	
	# check if the gene.loci
	if(!is(gene.loci, "GRanges")){
		stop("gene.loci is not genomic range object")
	}else{
		if((x=width(gene.loci)) == 0L){
			stop("length of gene flanking window is 0");
		}
	}
	
	idy = queryHits(findOverlaps(peak.use, gene.loci));
	if((x=length(idy))==0L){
		stop("no peaks are found within the gene flanking window")
	}else{
		data.use = data.use[,idy];
		peak.use = peak.use[idy];
		# remove peaks have low coverage
		idy = which(Matrix::colSums(data.use) > 0);
		if((x=length(idy)) == 0L){
			stop("no peaks are found within the gene flanking window")
		}
		data.use = data.use[,idy];
		peak.use = peak.use[idy];
	}
	
	
	cat("Epoch: performing statitical test ... \n", file = stderr())
	if(do.par){
	    # input checking for parallel options
		if(num.cores > 1){
	        if (num.cores == 1) {
	          num.cores = 1
	        } else if (num.cores > detectCores()) {
	          num.cores <- detectCores() - 1
	          warning(paste0("num.cores set greater than number of available cores(", parallel::detectCores(), "). Setting num.cores to ", num.cores, "."))
	        }
	      } else if (num.cores != 1) {
	        num.cores <- 1
		}
		
		cl <- makeCluster(num.cores, type = "SOCK");
	    registerDoSNOW(cl);
	    clusterEvalQ(cl, library(stats));
		peaks.id = seq(ncol(data.use));
		clusterExport(cl, c("data.use", "gene.val"), envir=environment());
		opts <- list(preschedule=TRUE);
		clusterSetRNGStream(cl, 10);
		models <- suppressWarnings(llply(.data=peaks.id, .fun=function(t) summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",], .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE));
		models <- do.call(rbind, models);
		stopCluster(cl);
	}else{
		peaks.id = seq(ncol(data.use));
		models = lapply(peaks.id, function(t){
			summary(stats::glm(y ~ x, data = data.frame(y=data.use[,t], x=gene.val), family = binomial(link='logit')))[["coefficients"]]["x",]
		})
		models <- do.call(rbind, models);
	}
	
	models[models[,"z value"] < 0,"Pr(>|z|)"] = 1;
	peak.use$beta = models[,"Estimate"];
	peak.use$zvalue = models[,"z value"];
	peak.use$stde = models[,"Std. Error"];
	peak.use$Pval = models[,"Pr(>|z|)"];
	peak.use$logPval = -log10(peak.use$Pval);
	cat("Epoch: Done ... \n", file = stderr())
	return(peak.use);
}
