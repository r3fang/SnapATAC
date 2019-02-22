#' Motif Analysis Using chromVAR
#' @export
runChromVAR <- function(object, ...) {
  UseMethod("runChromVAR", object);
}

#' @export
runChromVAR.default <- function(
	object,
	mat=c("pmat", "bmat"),
	genome=NULL,
	species="Homo sapiens",
	collection=c("CORE", "CNE", "PHYLOFACTS", "SPLICE", "POLII", "FAM", "PBM", "PBM_HOMEO", "PBM_HLH"),
	...
	){
		if(missing(object)){
			stop("'object' is missing");
		}
		
		if(class(object) != "snap"){
			stop("'object' is not a snap object");
		}
		
		
		mat = match.arg(mat);
		if(mat == "pmat"){			
			peaks.gr = object@peak;
			x = object@pmat;
		}else if(mat == "bmat"){
			peaks.gr = object@feature;			
			x = object@bmat;
		}else{
			stop("mat must be 'bmat' or 'pmat'")
		}
				
		rse <- SummarizedExperiment(assays = list(counts = t(x)), rowRanges = peaks.gr, colData = DataFrame(Cell_Type=1:nrow(x), depth=Matrix::rowSums(x)))
		rse <- addGCBias(rse, genome = genome);									  
		motifs <- getJasparMotifs(collection = collection, species=species, ...);
		motif_mm <- matchMotifs(motifs, rse, genome = genome);
		dev <- computeDeviations(object = rse, annotations = motif_mm);
		dev_mat = t(assay(dev));
		return(dev_mat);
	}