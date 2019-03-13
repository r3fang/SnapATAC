#' Call Peaks Using MACS2
#'
#' Identify peaks using selected cells. Fragments belonging to subset of cells 
#' are extracted and used to identify peaks using MACS2. This function requires
#' "MACS2" and "snaptools" preinstalled and excutable. 
#' 
#' @param obj A snap object.
#' @param file A snap file.
#' @param output.prefix Prefix of output file which will be used to generate output file names.
#' @param path.to.snaptools Path to snaptools excutable file.
#' @param path.to.macs Path to macs2 excutable file.
#' @param gsize effective genome size. 'hs' for human, 'mm' for mouse, 'ce' for C. elegans, 'dm' for fruitfly (default: None)
#' @param buffer.size Buffer size for incrementally increasing internal array size to store reads alignment information. In
#' most cases, you don't have to change this parameter. However, if there are very high coverage dataset that each barcode has
#' more than 10000 fragments, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read snap files) [1000].
#' @param macs.options String indicate options you would like passed to macs2. strongly do not recommand to change unless you know what you are doing. the default is '--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits'.
#' @param tmp.folder Directory to store temporary files. If not given, snaptools will automatically generate a temporary location to store temporary files. 
#' @param keep.minimal Keep minimal version of output [TRUE].
#' @return Return a data.frame object that contains the peak information
#'
#' @export
runMACS <- function(
	obj, 
	file,
	output.prefix,
	path.to.snaptools,
	path.to.macs,
	gsize,
	buffer.size=500,
	macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
	tmp.folder=NULL,
	keep.minimal=TRUE
){
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
		if((x=length(obj@barcode))==0L){
			stop("obj@barcode is empty");			
		}
		barcode.use = obj@barcode;
	}

	if(missing(file)){
		stop("file is missing");
	}else{
		if(!file.exists(file)){
			stop("file does not exist");
		}
		if(!isSnapFile(file)){
			stop("file is not a snap file");
		}
	}
	
	cat("Epoch: checking selected barcodes exist in snap file ... \n", file = stderr())
	barcode.flag = barcodeInSnapFile(barcode.use, file);
	if(any(barcode.flag == FALSE)){
		cat("the following barcode does not exist in snap file\n")
		cat(paste(barcode.use, "\n"));
		stop()
	}

	if(missing(output.prefix)){
		stop("output.prefix is missing");
	}
	
	if(missing(path.to.snaptools)){
		stop("path.to.snaptools is missing");
	}else{
		if(!file.exists(path.to.snaptools)){
			stop("path.to.snaptools does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.snaptools);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.snaptools is not an excutable file");
		}
		
	}

	if(missing(path.to.macs)){
		stop("path.to.macs is missing");
	}else{
		if(!file.exists(path.to.macs)){
			stop("path.to.macs does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.macs);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.macs is not an excutable file");
		}
		
	}
	
	if(missing(gsize)){
		stop("gsize is missing");
	}
	
	if(is.null(tmp.folder)){
		tmp.folder =  tempdir();
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist")
		}
	}
	
	# write down the barcode info
	tmp.file = tempfile(pattern = "run_macs_barcode", tmpdir = tmp.folder, fileext = ".txt")
	write.table(barcode.use, file = tmp.file, append = FALSE, quote = FALSE, sep = "\t",
	                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
	                 col.names = FALSE, qmethod = c("escape", "double"),
	                 fileEncoding = "")
	
	cat("Epoch: running macs2 for peak calling ...\n", file = stderr())
	flag = system2(command=path.to.snaptools, 
		args=c("call-peak", 
			   "--snap-file", file, 
			   "--barcode-file", tmp.file, 
			   "--output-prefix", output.prefix,
			   "--path-to-macs", dirname(path.to.macs),
			   "--gsize", gsize,
			   "--buffer-size", buffer.size,
			   "--macs-options", shQuote(macs.options),
			   "--tmp-folder", tmp.folder
			   )		
		);		

	if (flag != 0) {
	   	stop("'rumMACS' call failed");
	}				
	
	if(keep.minimal){
		system(paste("rm ", output.prefix, "_control_lambda.bdg", sep=""));
		system(paste("rm ", output.prefix, "_peaks.xls", sep=""));
		system(paste("rm ", output.prefix, "_summits.bed", sep=""));
	}
	return(read.table(paste(output.prefix, "_peaks.narrowPeak", sep="")));
}



#' Call Peaks Using MACS2 For All Clusters
#'
#' Identify peaks for all clusters. Fragments belonging to each subset or cluster of cells 
#' are extracted and used to identify peaks using MACS2. This function requires
#' "MACS2" and "snaptools" preinstalled and excutable. 
#' 
#' @param obj A snap object.
#' @param file A snap file.
#' @param output.prefix Prefix of output file which will be used to generate output file names.
#' @param path.to.snaptools Path to snaptools excutable file.
#' @param path.to.macs Path to macs2 excutable file.
#' @param min.cells min number of cells to perform peak calling [100]. Clusters with cells less 
#' than num.cells will be excluded.
#' @param num.cores number of cpus to use [1].
#' @param gsize effective genome size. 'hs' for human, 'mm' for mouse, 'ce' for C. elegans, 'dm' for fruitfly (default: None)
#' @param buffer.size Buffer size for incrementally increasing internal array size to store reads alignment information. In
#' most cases, you don't have to change this parameter. However, if there are very high coverage dataset that each barcode has
#' more than 10000 fragments, it's recommended to specify a smaller buffer size in order to decrease memory usage (but it will take longer time to read snap files) [1000].
#' @param macs.options String indicate options you would like passed to macs2. strongly do not recommand to change unless you know what you are doing. the default is '--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits'.
#' @param tmp.folder Directory to store temporary files. If not given, snaptools will automatically generate a temporary location to store temporary files. 
#' @param keep.minimal Keep minimal version of output [TRUE].
#' @return Return a GRanges object that contains the non-overlapping combined peaks
#' @importFrom parallel mclapply 
#' @importFrom GenomicRanges GRanges reduce
#' @export
runMACSForAll <- function(
	obj, 
	file,
	output.prefix,
	path.to.snaptools,
	path.to.macs,
	gsize,
	num.cores=1,
	min.cells=100,
	buffer.size=500,
	macs.options="--nomodel --shift 37 --ext 73 --qval 1e-2 -B --SPMR --call-summits",
	tmp.folder=NULL,
	keep.minimal=TRUE
){
	cat("Epoch: checking input parameters ... \n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is.snap(obj)){
			stop("obj is not a snap object");
		}		
		if((x=nrow(obj))==0L){
			stop("obj is empty");
		}
		if((x=length(obj@barcode))==0L){
			stop("obj@barcode is empty");			
		}
		barcode.use = obj@barcode;
		nclusters = length(levels(obj@cluster));
	}
	
	if(nclusters == 0L){
		stop("obj does not have cluster, runCluster first")
	}

	if(missing(file)){
		stop("file is missing");
	}else{
		if(!file.exists(file)){
			stop("file does not exist");
		}
		if(!isSnapFile(file)){
			stop("file is not a snap file");
		}
	}
	
	cat("Epoch: checking selected barcodes exist in snap file ... \n", file = stderr())
	barcode.flag = barcodeInSnapFile(barcode.use, file);
	if(any(barcode.flag == FALSE)){
		cat("the following barcode does not exist in snap file\n")
		cat(paste(barcode.use, "\n"));
		stop()
	}

	if(missing(output.prefix)){
		stop("output.prefix is missing");
	}
	
	if(missing(path.to.snaptools)){
		stop("path.to.snaptools is missing");
	}else{
		if(!file.exists(path.to.snaptools)){
			stop("path.to.snaptools does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.snaptools);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.snaptools is not an excutable file");
		}
		
	}

	if(missing(path.to.macs)){
		stop("path.to.macs is missing");
	}else{
		if(!file.exists(path.to.macs)){
			stop("path.to.macs does not exist");
		}
		
		flag = tryCatch({
			file_test('-x', path.to.macs);	
		},
		error=function(cond){
			return(FALSE)
		})
		if(flag == FALSE){
			stop("path.to.macs is not an excutable file");
		}
		
	}
	
	if(missing(gsize)){
		stop("gsize is missing");
	}
	
	if(is.null(tmp.folder)){
		tmp.folder =  tempdir();
	}else{
		if(!dir.exists(tmp.folder)){
			stop("tmp.folder does not exist")
		}
	}
	
    peak.ls <-  parallel::mclapply(as.list(levels(obj@cluster)), function(x){
       num.cells = length(which(obj@cluster == x));
       if(num.cells < min.cells){
           return(GenomicRanges::GRanges())
       }
       peaks.df = runMACS(
       obj=obj[which(obj@cluster==x),], 
       file=file, 
       output.prefix=paste(output.prefix, x, sep="."),
       path.to.snaptools=path.to.snaptools,
       path.to.macs=path.to.macs,
       gsize=gsize, 
       buffer.size=buffer.size, 
       macs.options=macs.options,
       tmp.folder=tmp.folder
       );
       peaks.gr = GenomicRanges::GRanges(peaks.df[,1], IRanges(peaks.df[,2], peaks.df[,3]))
       peaks.gr
     }, mc.cores=num.cores)
	
	peaks.gr = GenomicRanges::reduce(do.call(c, peak.ls));
	return(peaks.gr);
}


