#' Estimate Motif Variability Using ChromVAR
#'
#' @param obj A snap object.
#' @param input.mat Input matrix used for chromVAR analysis c("pmat", "bmat").
#' @param genome BSgenome object that contains the sequence for the corresponding genome (i.e. BSgenome.Hsapiens.UCSC.hg19).
#' @param min.count Min count for features. Features of count lower than min.count will be filtered [10].
#' @param species Species or jaspar code used by 'getJasparMotifs'to get the motif database ["Homo sapiens"].
#' @return Return a matrix object that contains the cell-by-motif matrix
#'
#' @export
runChromVAR <- function(
	obj, 
	input.mat=c("pmat", "bmat"),
	genome=BSgenome.Hsapiens.UCSC.hg19,
	min.count=10,
	species="Homo sapiens"
){
	cat("Epoch: checking depedent packages ... \n", file = stderr())
	if (!(requireNamespace("SummarizedExperiment", quietly = TRUE))) {
	      stop("Please install package 'SummarizedExperiment'")
	  }

  	if (!(requireNamespace("chromVAR", quietly = TRUE))) {
  	      stop("Please install package 'chromVAR'");
	}
	
	if (!(requireNamespace("motifmatchr", quietly = TRUE))) {
		stop("Please install package 'motifmatchr'")
	}
	
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}		
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
	if(!is(min.count, "numeric")){
		stop("min.count is not a numeric");
	}else{
		min.count = max(min.count, 0);
	}
	
	if(!(is(genome, "BSgenome"))){
		stop("genome is not a BSgenome object");
	}
	
	if(!(is(species, "character"))){
		stop("species is not a character")
	}
	
	idy = which(Matrix::colSums(data.use) >= min.count);
	data.use = data.use[,idy,dropping=TRUE]
	
	if(dim(data.use)[2] == 0){
		stop("input matrix is empty after filering low-coverage features")
	}
	peak.use = peak.use[idy];
	
	cat("Epoch: creating chromVAR object ... \n", file = stderr())
	rse <- SummarizedExperiment(
		assays = list(counts = t(data.use)), 
				 rowRanges = peak.use, 
				 colData = DataFrame(Cell_Type=1:nrow(data.use), depth=Matrix::rowSums(data.use))
	);
	cat("Epoch: computing GC bias ... \n", file = stderr());
	rse <- addGCBias(rse, genome = genome);
	cat("Epoch: getting JASPAR motifs ... \n", file = stderr());
	motifs <- getJasparMotifs(collection = "CORE", species=species);
	cat("Epoch: scanning motifs in the peaks ... \n", file = stderr());
	motif_mm <- matchMotifs(motifs, rse, genome = genome);
	cat("Epoch: calculating motif variability between cells ... \n", file = stderr());
	dev <- computeDeviations(object = rse, annotations = motif_mm);
	dev_mat = t(assay(dev));
	cat("Epoch: Done ... \n", file = stderr());
	return(dev_mat);
}
