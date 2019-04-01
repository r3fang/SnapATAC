#' Motif Analysis Using Homer
#' 
#' Program will find de novo and known motifs in regions in the genome using HOMER
#' 
#' @param obj A snap object.
#' @param result.dir Directory to store Homer results.
#' @param mat matrix to use c("pmat", "bmat").
#' @param path.to.homer Directory path to "findMotifsGenome.pl" excutable file.
#' @param genome Genome hg19 for human and mm10 for mouse.
#' @param num.cores Number of cores to use [10].
#' @param motif.length Motif length (default=8,10,12). NOTE: values greater 12 may cause the program to run out of memory.
#' @param scan.size fragment size to use for motif finding [200].
#' @param optimize.count global optimization: searches for strings with # mismatches [2].
#' @param background Genomic positions to be used as background. Removes background positions overlapping with target positions [automatic]
#' @param local.background Wehther to use local background [FALSE]
#' @param only.known Only to search for known motifs [TRUE]
#' @param only.denovo Only to search for de novo motifs [FALSE].
#' @param fdr.num Calculate empirical FDR for de novo discovery #=number of randomizations [5].
#' @param cache Size in MB for statistics cache [500].
#' @param keep.minimal Keep minimal version of output [FALSE].
#' @param ... Arguments passed to "findMotifsGenome.pl".
#' @importFrom utils file_test write.table read.csv 
#' @export
runHomer<- function(
	obj,
	result.dir,
	mat,
	path.to.homer,
	genome,
	num.cores,
	motif.length,
	scan.size,
	optimize.count,
	background,
	local.background,
	only.known,
	only.denovo,
	fdr.num,
	cache,
	keep.minimal,
	...
	) {
  UseMethod("runHomer", obj);
}

#' @export
runHomer.default <- function(
	obj,
	result.dir,
	mat=c("pmat", "bmat"),
	path.to.homer=NULL,
	genome = 'mm10',
	num.cores = 10,
	motif.length = 10,
	scan.size = 200,
	optimize.count = 2,
	background = 'automatic',
	local.background = FALSE,
	only.known = TRUE,
	only.denovo = FALSE,
	fdr.num = 5,
	cache = 100,
	keep.minimal = FALSE,
	...
	){
		path.to.homer <- normalizePath(path.to.homer);
		if (!file_test('-x', path.to.homer)) {
			stop(path.to.homer, " does not exist or is not executable; check your path.to.homer parameter")
		}		
		
		if(missing(result.dir)){
			stop("result.dir is missing")
		}else{
			if(dir.exists(result.dir)){
				stop("result.dir already exists, remove it first");			
			}
		}

		# check input
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}
		
		mat = match.arg(mat);
		if(mat == "pmat"){
			x <- as.data.frame(obj@peak)[,c(1:3)];
		}else if(mat == "bmat"){
			x <- as.data.frame(obj@feature)[,c(1:3)];
		}
		
		if(nrow(x) == 0){
			stop("feature for motif analysis is empty")
		}
		
	    ## Error checking -----------------------------------------------------
	    if (background != "automatic" && local.background != FALSE) {
	        stop("`background` and `local.background` are mutually exclusive; use only one")
	    }
		
	    if (only.known != FALSE & only.denovo != FALSE) {
	        stop("Both `only.known` and `only.denovo` set to `TRUE`; pick one")
	    }
        
		
		path.to.homer <- normalizePath(path.to.homer);
		if (!file_test('-x', path.to.homer)) {
			stop(path.to.homer, " does not exist or is not executable; check your path.to.homer parameter")
		}		
		
		
	    if ("data.frame" %in% class(x)) {
	        target_bed <- tempfile("target_", tmpdir=result.dir)
			write.table(x, file=target_bed, row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
	    } else {
	        if (file.exists(x) != TRUE) {
	            stop("Check that your bed file for `x` exists")
	        }        
	        target_bed <- x
	    }
		
	    if (!("automatic" %in% background)) {
	        if ("data.frame" %in% class(background)) {
	            background_bed <- tempfile("background_", tmpdir=result.dir)
	            #.write_bed(background, path = background_bed)
				write.table(target_bed, file=background_bed, row.names=FALSE, col.names=FALSE, sep="\t", quote = FALSE)
	        } else {
	            if (file.exists(background) != TRUE) {
	                stop("Check that your bed file for `background` exists")
	            }        
	            background_bed <- background
	        }
	    }

	    ## Make HOMER results output dir
	    system(paste("mkdir -p", result.dir))
        
	    cmd <- paste(
	        path.to.homer,
	        target_bed, genome, result.dir,
	        "-len", paste0(motif.length, collapse = ","),
	        "-size", scan.size,
	        "-S", optimize.count,
	        "-p", num.cores,
	        "-cache", cache,
	        "-fdr", fdr.num
	    )
		
	    if (!("automatic" %in% background)) {
	        cmd <- paste(cmd, "-bg", background_bed)
	    }
	    if (local.background != FALSE) {
	        cmd <- paste(cmd, "-local", local.background)
	    }
	    if (only.known == TRUE) {
	        cmd <- paste(cmd, "-nomotif")
	    }
	    if (only.denovo == TRUE) {
	        cmd <- paste(cmd, "-noknown")
	    }
	    if (scan.size == "given") {
	        cmd <- paste(cmd, "-chopify")
	    }
		
		system(cmd);
		
	    ## Remove extraneous files if desired
	    if (keep.minimal == TRUE) {
	        extra_files <- c("homerResults.html",
	                         "knownResults.html",
	                         "homerMotifs.motifs*",
	                         "motifFindingParameters.txt",
	                         "seq.autonorm.tsv",
	                         "*tmp*")
	        extra_dirs <- c("homerResults",
	                        "knownResults",
	                        "randomizations")
	        remove_extra <- paste(c(paste0("rm -f ", result.dir, "/", extra_files),
	                                paste0("rm -Rf ", result.dir, "/", extra_dirs)),
	                              collapse = "; ")
	        system(remove_extra)
	    }
		x = read.csv(paste0(result.dir, "/knownResults.txt"), sep="\t", header=TRUE);
	    system("rm -f *.tmp");
		return(x)
}

	

