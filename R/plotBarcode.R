#' Plot Snap Barocde Stats
#'
#' Plot the basic barcode statistics
#'
#' @param obj a snap object
#' @examples
#' snap_file <- system.file("extdata", "Fang.3C1.snap", package = "nucleus");
#' showBinSizes(snap_file);
#' x.sp <- readSnap(snap_file, binSize=5000, metaData=TRUE);
#' plotBarcode(x.sp);
#' @export
plotBarcode <- function(obj) {
  UseMethod("plotBarcode", obj);
}

#' @export
plotBarcode.default <- function(obj){	
	# check the input
	if(!(class(obj)=="snap")){stop(paste("Error @plotBarcode: obj is not a snap object!", sep=""))};
	if(nrow(obj@metaData) == 0){stop(paste("Error @plotBarcode: obj@metaData is empty!", sep=""))}
	barcode = getMetaData(obj);
	par(mfrow=c(3,2));
	hist(log(barcode$TN + 1,10), col="grey", main="log10(Total Fragments)", xlab="log10(Total Fragments + 1)", breaks=50);
	hist(log(barcode$UQ + 1,10), col="grey", main="log10(UMI)", xlab="log10(UMI + 1)", breaks=50);
	hist((barcode$UM+1)/(barcode$TN+1), col="grey", main="Mappability Ratio", xlab="Mappability", breaks=50, xlim=c(0, 1));
	hist((barcode$PP+1)/(barcode$UM+1), col="grey", main="Proper Paired Ratio", xlab="Proper Paired", breaks=50, xlim=c(0, 1));
	hist(1 - (barcode$UQ+1)/(barcode$PP+1), col="grey", main="Duplicate Rate Ratio", xlab="Duplicate Rate", breaks=50, xlim=c(0, 1));
	hist((barcode$CM+1) / (barcode$UQ+1), col="grey", main="chrM Rate Ratio", xlab="chrM Rate", breaks=50, xlim=c(0, 1));
}

