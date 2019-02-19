#' Combine snap objects
#'
#' Takes two snap objects and combines them.
#'
#' @param obj1 a snap object
#' @param obj2 a snap object
#' @return a combined snap object
#' @examples
#' file <- system.file("extdata", "Fang.3C1.snap", package = "SNAPATAC");
#' x1.sp <- readSnap(file, binSize=5000, metaData=TRUE);
#' file <- system.file("extdata", "Fang.3C2.snap", package = "SNAPATAC");
#' x2.sp <- readSnap(file, binSize=5000, metaData=TRUE);
#' rBind(x1.sp, x2.sp);
#' @export
rBind <- function(obj1, ...){
  UseMethod("rBind", obj1);
}

rBind.default <- function(obj1, obj2){
	# only the following slots can be combined
	# barcode, feature, metaData, cmat, bmat
	# among these slots, barcode, feature, cmat are enforced, the others are optional
	# the rest slots must be set to be empty
	
	if(!is.snap(obj1)){stop(paste("Error @rBind: obj1 is not a snap object!", sep=""))};
	if(!is.snap(obj2)){stop(paste("Error @rBind: obj2 is not a snap object!", sep=""))};

	# barcode from obj1 and obj2
	barcode1 = obj1@barcode;
	barcode2 = obj2@barcode;	
	
	# check barcode name, if there exists duplicate barcode raise error and exist
	if(length(unique(c(barcode1, barcode2))) < length(barcode1) + length(barcode2)){
		stop("Error: @rBind: identifcal barcodes found in obj1 and obj2!")
	}
	barcode = c(barcode1, barcode2);
	
	# check meta data
	if(nrow(obj1@metaData) > 0 && nrow(obj2@metaData) > 0){
		metaData = rbind(obj1@metaData, obj2@metaData);		
	}else{
		metaData = data.frame();
	}
	
	# check feature
	feature1 = obj1@feature;
	feature2 = obj2@feature;
	if((length(feature1) == 0) != (length(feature2) == 0)){
		stop("different feature found in obj1 and obj2!")
	}else{
		if(length(feature1) > 0){
			if(FALSE %in% (feature1$name == feature2$name)){
				stop("Error: @rBind: different feature found in obj1 and obj2!")
			}
			feature = feature1;					
		}else{
			feature = feature1;								
		}
	}
	
	# check peak
	peak1 = obj1@peak;
	peak2 = obj2@peak;
	if((length(peak1) == 0) != (length(peak2) == 0)){
		stop("different peak found in obj1 and obj2!")
	}else{
		if(length(peak1) > 0){
			if(FALSE %in% (peak1$name == peak2$name)){
				stop("Error: @rBind: different feature found in obj1 and obj2!")
			}
			peak = peak1;					
		}else{
			peak = peak1;								
		}
	}
	
	# check bmat	
	bmat1 = obj1@bmat;
	bmat2 = obj2@bmat;
	if((length(bmat1) == 0) != (length(bmat2) == 0)){
		stop("bmat has different dimentions in obj1 and obj2!")
	}else{
		bmat = Matrix::rBind(bmat1, bmat2);
	}

	# check gmat	
	gmat1 = obj1@gmat;
	gmat2 = obj2@gmat;
	if((length(gmat1) == 0) != (length(gmat2) == 0)){
		stop("gmat has different dimentions in obj1 and obj2!")
	}else{
		gmat = Matrix::rBind(gmat1, gmat2);
	}

	# check pmat	
	pmat1 = obj1@pmat;
	pmat2 = obj2@pmat;
	if((length(pmat1) == 0) != (length(pmat2) == 0)){
		stop("pmat has different dimentions in obj1 and obj2!")
	}else{
		pmat = Matrix::rBind(pmat1, pmat2);
	}

	res = newSnap();
	res@barcode = barcode;
	res@metaData = metaData;
	res@bmat = bmat;
	res@pmat = pmat;
	res@feature = feature;
	res@peak = peak;
	res@gmat = gmat;
	return(res)
}






