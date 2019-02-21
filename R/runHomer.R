#' Motif Analysis Using Homer
#' @export

runHomer<- function(object, ...) {
  UseMethod("runHomer", object);
}

#' @export
runHomer.default <- function(
	object,
	mat=c("pmat", "bmat"),
	result_dir=NULL,
	genome = 'mm10',
	cores = 10,
	motif_length = 10,
	scan_size = 300,
	optimize_count = 2,
	background = 'automatic',
	local_background = FALSE,
	only_known = FALSE,
	only_denovo = FALSE,
	fdr_num = 5,
	cache = 100,
	overwrite = TRUE,
	keep_minimal = FALSE,
	...
	){
		if (is.null(result_dir)) {
			result_dir <- tempfile(pattern='homer')
		}
		
		# check input
		if(class(object) != "snap"){
			stop("'object' is not a snap object")
		}
		
		mat = match.arg(mat);
		if(mat == "pmat"){
			peaks <- as_tibble(as.data.frame(object@peak)[,c(1:3)]);
		}else if(mat == "bmat"){
			peaks <- as_tibble(as.data.frame(object@feature)[,c(1:3)]);
		}
		
		if(nrow(peaks) == 0){
			stop("feature for motif analysis is empty")
		}
		
	    ## Error checking -----------------------------------------------------
	    if (overwrite == FALSE & dir.exists(path)) {
	        stop("Output directory exists (set `overwrite = TRUE` to bypass)")
	    }

	    if (background != "automatic" && local_background != FALSE) {
	        stop("`background` and `local_background` are mutually exclusive; use only one")
	    }
	    if (only_known != FALSE & only_denovo != FALSE) {
	        stop("Both `only_known` and `only_denovo` set to `TRUE`; pick one")
	    }
        
	    ## File format for `x` and `background` (bed vs data.frame)
	    ## Write bed files if data.frame; assign if a character
	    if ("data.frame" %in% class(x)) {
	        target_bed <- tempfile("target_")
			write.table(x, file=target_bed, row.name=FALSE, col.name=FALSE, sep="\t", quote = FALSE)
	    } else {
	        if (file.exists(x) != TRUE) {
	            stop("Check that your bed file for `x` exists")
	        }        
	        target_bed <- x
	    }
		
	    if (!("automatic" %in% background)) {
	        if ("data.frame" %in% class(background)) {
	            background_bed <- tempfile("background_")
	            #.write_bed(background, path = background_bed)
				write.table(target_bed, file=background_bed, row.name=FALSE, col.name=FALSE, sep="\t", quote = FALSE)
	        } else {
	            if (file.exists(background) != TRUE) {
	                stop("Check that your bed file for `background` exists")
	            }        
	            background_bed <- background
	        }
	    }

		
	    ## Make HOMER results output dir
	    system(paste("mkdir -p", result_dir))
    
	    ## Construct and run command
	    homer_base <- get_homer_base()
    
	    cmd <- paste(
	        "/mnt/silencer2/home/r3fang/apps/homer/bin/findMotifsGenome.pl",
	        target_bed, genome, result_dir,
	        "-len", paste0(motif_length, collapse = ","),
	        "-size", scan_size,
	        "-S", optimize_count,
	        "-p", cores,
	        "-cache", cache,
	        "-fdr", fdr_num
	    )
		
	    if (!("automatic" %in% background)) {
	        cmd <- paste(cmd, "-bg", background_bed)
	    }
	    if (local_background != FALSE) {
	        cmd <- paste(cmd, "-local", local_background)
	    }
	    if (only_known == TRUE) {
	        cmd <- paste(cmd, "-nomotif")
	    }
	    if (only_denovo == TRUE) {
	        cmd <- paste(cmd, "-noknown")
	    }
	    if (scan_size == "given") {
	        cmd <- paste(cmd, "-chopify")
	    }
	    system(cmd)

	    ## Remove extraneous files if desired
	    if (keep_minimal == TRUE) {
	        extra_files <- c("homerResults.html",
	                         "knownResults.html",
	                         "homerMotifs.motifs*",
	                         "motifFindingParameters.txt",
	                         "seq.autonorm.tsv",
	                         "*tmp*")
	        extra_dirs <- c("homerResults",
	                        "knownResults",
	                        "randomizations")
	        remove_extra <- paste(c(paste0("rm -f ", result_dir, "/", extra_files),
	                                paste0("rm -Rf ", result_dir, "/", extra_dirs)),
	                              collapse = "; ")
	        system(remove_extra)
	    }
		
		x = read.csv(paste0(result_dir, "/knownResults.txt"), sep="\t", head=TRUE);
	    system("rm -f *.tmp");
		return(x)
}

	

