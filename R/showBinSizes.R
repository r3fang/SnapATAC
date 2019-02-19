#' Snap-format file
#'
#' This function takes a snap-format file name as input and check the bin 
#' sizes or resolutions have been generated for count matrix
#'
#' @param file character. input snap-format file name
#'
#' @return integer vector. A vector of integers indicating the bin sizes
#'
#' @examples
#' fname <- system.file("extdata", "Schep.snap", package = "nucleus")
#' showBinSizes(fname)
#'
#' @export
#'
showBinSizes <- function(file) {
  UseMethod("showBinSizes", file);
}

#' @export
showBinSizes.default <- function(file){
	# this is to check what are the binSizes have been generated
	if(!file.exists(file)){stop("file does not exist!")}
	if(!isSnapFile(file)){stop("input file is not a snap file\n")}		
	binSizeList = tryCatch(binSizeList <- h5read(file, '/AM/binSizeList'), error = function(e) {print(paste("Warning @check.bin.size: '/AM/binSizeList' not found in ", file)); return(numeric())})
	return(binSizeList);
}
