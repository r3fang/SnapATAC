#' Read meta data from a snap file
#'
#' Take a snap file as input and read the barcode session only.
#'
#' @param file character for the snap-format file name which the data are to be read from.
#' @return A data frame contains barcodes and its attributes
#'
#' @examples
#' snap_file <- system.file("extdata", "Fang.3C1.snap", package = "SNAPATAC");
#' barcodes.df <- readMetaData(snap_file);
#' @export
readMetaData <- function(file, ...) {
  UseMethod("readMetaData", file);
}

#' @export
readMetaData.default <- function(file){
	if(!file.exists(file)){stop(paste("Error @readMetaData: ", file, " does not exist!", sep=""))};
	if(!isSnapFile(file)){stop(paste("Error @readMetaData: ", file, " is not a snap-format file!", sep=""))};
	# close the previously opened H5 file
	H5close();
	barcode = as.character(tryCatch(barcode <- h5read(file, '/BD/name'), error = function(e) {print(paste("Warning @readSnap: 'BD/name' not found in ",file)); return(vector(mode="character", length=0))}));
	TN = as.numeric(tryCatch(TN <- h5read(file, '/BD/TN'), error = function(e) {print(paste("Warning @readMetaData: 'BD/TN' not found in ",file)); return(c())}));
	UM = as.numeric(tryCatch(UM <- h5read(file, '/BD/UM'), error = function(e) {print(paste("Warning @readMetaData: 'BD/UM' not found in ",file)); return(c())}));
	PP = as.numeric(tryCatch(PP <- h5read(file, '/BD/PP'), error = function(e) {print(paste("Warning @readMetaData: 'BD/PP' not found in ",file)); return(c())}));
	UQ = as.numeric(tryCatch(UQ <- h5read(file, '/BD/UQ'), error = function(e) {print(paste("Warning @readMetaData: 'BD/UQ' not found in ",file)); return(c())}));
	CM = as.numeric(tryCatch(CM <- h5read(file, '/BD/CM'), error = function(e) {print(paste("Warning @readMetaData: 'BD/CM' not found in ",file)); return(c())}));
	data.frame(barcode, TN, UM, PP, UQ, CM)
}
