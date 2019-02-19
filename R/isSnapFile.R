#' Snap
#' 
#' This function takes a file name as input and check if the file
#' is a snap-formated file 
#'
#' @param file character; A character indicates the file name
#' 
#' @return Weather a given file is a snap-formatted file
#'
#' @examples
#' fname <- system.file("extdata", "Schep.snap", package = "nucleus")
#' isSnapFile(fname)
#' 
#' @export
#' 
isSnapFile <- function(file) {
  UseMethod("isSnapFile");
}

#' @describeIn isSnapFile Default Interface
#' @export
isSnapFile.default <- function(file){
	# this is to check if a file is a snap file;

	# return TRUE if fname is a mona file, otherwise FALSE
	if(!file.exists(file)){ stop("Error @isSnapFile: does not exist!")}

	# check if the file is a H5 file
	monals = tryCatch(
		monals <- h5ls(file), 
		error = function(e) {
			return(character())
	})
	
	if(length(monals) != 0){
		magicString = as.character(tryCatch(magicString <- h5read(file, '/HD/MG'), error = function(e) {print(paste("Warning @isSnapFile: 'HD/MG' not found in ", fname)); return(character())}))
		if(length(magicString)==0){
			return(FALSE);
		}
		if(magicString == "SNAP" || magicString == "MONA"){
			return(TRUE);
		}
	}
	return(FALSE);
}

