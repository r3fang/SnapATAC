#' Plot Barocde Quality Control Distribution
#'
#' @param obj a snap object
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @param ... Arguments passed to hist method.
#' @examples 
#' data(demo.sp);
#' plotBarcode(demo.sp, col="grey", border="grey");
#' @importFrom graphics par hist 
#' @importFrom grDevices pdf dev.off
#' @export
plotBarcode <- function(obj, pdf.file.name, pdf.width, pdf.height, ...) {
  UseMethod("plotBarcode", obj);
}

#' @export
plotBarcode.default <- function(
	obj,
	pdf.file.name=NULL,
	pdf.width=7,
	pdf.height=7,
	...
	){	

	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object");
		}
		ncell = nrow(obj);
		if(nrow(obj@metaData) == 0){stop(paste("obj deos not have meta data", sep=""))}
		# check the input
		barcode = obj@metaData;
	}
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
			
	par(mfrow=c(3,2));
	hist(log(barcode$TN + 1,10), main="log10(Total Fragments)", xlab="log10(Total Fragments + 1)", ...);
	hist(log(barcode$UQ + 1,10), main="log10(UMI)", xlab="log10(UMI + 1)", ...);
	hist((barcode$UM+1)/(barcode$TN+1), main="Mappability Ratio", xlab="Mappability", xlim=c(0, 1), ...);
	hist((barcode$PP+1)/(barcode$UM+1), main="Proper Paired Ratio", xlab="Proper Paired",xlim=c(0, 1), ...);
	hist(1 - (barcode$UQ+1)/(barcode$PP+1), main="Duplicate Rate Ratio", xlab="Duplicate Rate", xlim=c(0, 1), ...);
	hist((barcode$CM+1) / (barcode$UQ+1), main="chrM Rate Ratio", xlab="chrM Rate", xlim=c(0, 1), ...);	
	if(!is.null(pdf.file.name)){
		dev.off()		
	}
	graphics::par(mfrow=c(1,1));
}

#' Plot Bin Coverage Distribution
#'
#' @param obj a snap object
#' @param rm.zeros Remove bins of zero coverage when ploting the coverage distribution [TRUE].
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @param ... Arguments passed to hist method.
#' @examples 
#' data(demo.sp);
#' plotBinCoverage(demo.sp, col="grey", border="grey");
#' @importFrom stats sd
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics hist
#' @export
plotBinCoverage <- function(obj, rm.zeros, pdf.file.name, pdf.width, pdf.height, ...) {
  UseMethod("plotBinCoverage", obj);
}

#' @export
plotBinCoverage.default <- function(
	obj, 
	rm.zeros=TRUE,
	pdf.file.name=NULL,
	pdf.width=7, 
	pdf.height=7, 
	...
){	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object");
		}
		ncell = nrow(obj);
		data.use = obj@bmat;
		
		if((x=nrow(data.use)) == 0L){
			stop("cell-by-bin matrix is empty, run addBmatToSnap first")
		}
	}
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		grDevices::pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	
	cov = Matrix::colSums(data.use);
	idy = seq(cov);	
	cov = cov[which(cov > 0)];	
	cov = log10(cov + 1);
	cov = (cov - mean(cov)) / sd(cov);
	hist(cov, ...);
	
	if(!is.null(pdf.file.name)){
		grDevices::dev.off()		
	}
	graphics::par(mfrow=c(1,1));
}

#' Visulization
#'
#' @param obj A snap object.
#' @param method Visulization method c("tsne", "umap").
#' @param point.size Point size [1].
#' @param point.shape Point shape type [19].
#' @param point.alpha Point transparancy level [0.8].
#' @param point.color Color of point. Two options, points color by 
#' cell cluster label or sample ID c("cluster", "sample"). 
#' @param text.add Whether to add cluster label text at the centroid of each cluster [TRUE].
#' @param text.size Cluster label text size [5].
#' @param text.color Cluster label text color ["black"].
#' @param text.halo.add Add halo to cluster label text [TRUE].
#' @param text.halo.color Halo color ["white"].
#' @param text.halo.width Halo width [0.5].
#' @param legend.add Add a legend to the plot [FALSE].
#' @param legend.pos Position of the legend c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center").
#' @param legend.text.size Size of the text in the legend [1].
#' @param legend.text.color Color of the text in the legend ["black"].
#' @param down.sample Downsample the original cells to down.sample cells to ovoid large dataset [10,000].
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @param ... Arguments passed to plot method.
#'
#' @examples 
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	mat="bmat"
#'	);
#' demo.sp = runNormJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir()
#'	);
#' demo.sp = runDimReduct(
#'	obj=demo.sp, 
#'	pc.num=10, 
#'	input.mat="jmat"
#'	);
#' demo.sp = runKNN(
#'	obj=demo.sp, 
#'	pca.dims=1:5, 
#'	k=15, 
#'	snn=TRUE, 
#'	save.knn=FALSE
#'	);
#' demo.sp = runCluster(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	louvain.lib="R-igraph"
#'	);
#' demo.sp = runViz(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	pca.dims=1:5, 
#'	method="Rtsne"
#'	);
#' plotViz(
#' 	obj=demo.sp, 
#' 	method="tsne", 
#' 	point.color="sample", 
#' 	text.add=FALSE,
#'	legend.add=TRUE
#' );
#' 
#' @importFrom grDevices pdf dev.off
#' @importFrom methods slot
#' @importFrom scales alpha
#' @importFrom graphics plot text title legend
#' @export
plotViz <- function(obj, 
	method, 
	point.size, 
	point.shape, 
	point.alpha, 
	point.color, 
	text.add, 
	text.size, 
	text.color, 
	text.halo.add, 
	text.halo.color, 
	text.halo.width, 
	legend.add, 
	legend.pos, 
	legend.text.size,
	legend.text.color,
	down.sample, 
	pdf.file.name, 
	pdf.width, 
	pdf.height, 
	...
){
  UseMethod("plotViz", obj);
}

#' @export
plotViz.default <- function(obj, 
		method=c("tsne", "umap"), 
		point.size=1, 
		point.shape=19, 
		point.alpha=0.8, 
		point.color=c("cluster", "sample"),
		text.add=TRUE,
		text.size=1, 
		text.color="black",
		text.halo.add=TRUE,
		text.halo.color="white",
		text.halo.width=0.2,
		legend.add=FALSE,
		legend.pos=c("bottomleft", "bottom", "left", "topleft", "top", "topright", "right", "center"),
		legend.text.size=1,
		legend.text.color="black",
		down.sample=10000,
		pdf.file.name=NULL,
		pdf.width=7, 
		pdf.height=7,
		...
){	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object");
		}
		ncell = nrow(obj);
	}
		
	if(is.integer(down.sample)){
		stop("down.sample must be an integer");
	}
	
	method = match.arg(method);
	data.use = as.data.frame(slot(obj, method));
	if((x=nrow(data.use)) == 0L){
		stop("visulization method does not exist, run runViz first!")
	}

	cluster = slot(obj, point.color);
	if(length(cluster) == 0L){
		warning("cluster does not exist, text.add is ignored")
		text.add = FALSE;
	}
	point.color = match.arg(point.color);
	if(length(cluster) != 0L){	
		data.use$col = factor(cluster);
	}else{
		data.use$col = factor(1);
	}
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	legend.pos = match.arg(legend.pos);
	down.sample = min(down.sample, ncell);
	idx.ds = sort(sample(seq(ncell), down.sample));
	data.use = data.use[idx.ds,,drop=FALSE]
										
	colPanel = createColorPanel(length(unique(data.use$col)));
	graphics::plot(x=data.use[,1],
				   y=data.use[,2],
		 		   cex=point.size, 
		 		   pch=point.shape, 
		 		   col=scales::alpha(colPanel[factor(data.use$col)], point.alpha),
		 		   xlab="",
		 		   ylab="",
		 		   yaxt='n',
		 		   xaxt='n',
		 		   axes=FALSE,
		 		   ...
				   );
	graphics::title(ylab="Dim-2", line=0.5, cex.lab=1.2, font.lab=2)
	graphics::title(xlab="Dim-1", line=0.5, cex.lab=1.2, font.lab=2)
	graphics::box(lwd=2);		
  	
	if(text.add){
		xx = findCentrod(data.use[,c(1,2)], data.use$col);
		textHalo(x=xx[,1], y=xx[,2], labels = xx[,3], col=text.color, bg=text.halo.color, r=text.halo.width, cex=text.size);
  	}
	
	if(legend.add){
		legend(legend.pos, 
		  legend = levels(factor(data.use$col)), 
		  col = colPanel,
		  pch = point.shape,
		  pt.cex=1, 
		  bty = "n", 
		  cex = legend.text.size, 
		  text.col = legend.text.color,
		  horiz = FALSE
		  )		
	}

	if(!is.null(pdf.file.name)){
		dev.off()		
	}
}


#' Elbow Plot for Dimentionality Reduction Result
#'
#' @param obj A snap object
#' @param point.size Point size [2].
#' @param point.shape Point shape type [19].
#' @param point.color Point color ["red"].
#' @param point.alpha Point transparancy level [1].
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @param ... Arguments passed to plot function. 
#'
#' @examples 
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	mat="bmat"
#'	);
#' demo.sp = runNormJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir()
#'	);
#' demo.sp = runDimReduct(
#'	obj=demo.sp, 
#'	pc.num=10, 
#'	input.mat="jmat"
#'	);
#' plotDimReductElbow(demo.sp);
#' 
#' @importFrom grDevices pdf dev.off
#' @importFrom methods slot
#' @importFrom scales alpha
#' @importFrom graphics plot text title
#' @export
plotDimReductElbow <- function(obj, point.size, point.shape, point.color, point.alpha, pdf.file.name, pdf.height, pdf.width, ...){
  UseMethod("plotDimReductElbow", obj);
}

#' @export
plotDimReductElbow.default <- function(
	obj, 
	point.size=1.5,
	point.shape=19,
	point.color="red",
	point.alpha=1,
	pdf.file.name=NULL,
	pdf.height=7,
	pdf.width=7,
	...
){
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object");
		}
		ncell = nrow(obj);
		if(!isDimReductComplete(obj@smat)){
			stop("obj does not have valid dim.reduct object, run 'runDimReduct' first");			
		}
	}

	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	data.use = data.frame(PC=1:length(obj@smat@sdev), sd=obj@smat@sdev);	
	plot(x=data.use[,1], y=data.use[,2], cex=point.size, pch=point.shape, col=alpha(point.color, point.alpha), xlab="PCs", ylab="Standard Deviation of PCs");
		
	if(!is.null(pdf.file.name)){
		dev.off()
	}
}


#' Pairwise plot for Dimentionality Reduction Result
#'
#' @param obj A snap object
#' @param pca.dims PC dimetions to plot [1:30]
#' @param point.size Point size [2].
#' @param point.shape Point shape type [19].
#' @param point.color Point color ["grey"].
#' @param point.alpha Point transparancy level [1].
#' @param down.sample Number of cells to plot. down.sample cells will be randomly downsampled for plot if there are more cells.
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' 
#' @examples 
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	mat="bmat"
#'	);
#' demo.sp = runNormJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir()
#'	);
#' demo.sp = runDimReduct(
#'	obj=demo.sp, 
#'	pc.num=10, 
#'	input.mat="jmat"
#'	);
#' plotDimReductPW(demo.sp, pca.dims=1:10);
#' 
#' @importFrom grDevices pdf dev.off
#' @importFrom methods slot
#' @importFrom scales alpha
#' @importFrom graphics plot text title mtext
#' @export
plotDimReductPW <- function(obj, pca.dims, point.size, point.color, point.shape, point.alpha, down.sample, pdf.file.name, pdf.height, pdf.width){
  UseMethod("plotDimReductPW", obj);
}

#' @export
plotDimReductPW.default <- function(
	obj, 
	pca.dims=1:50,
	point.size=0.5,
	point.color="grey",
	point.shape=19,
	point.alpha=0.5,
	down.sample=3000,
	pdf.file.name=NULL, 
	pdf.height=7, 
	pdf.width=7
){
	
	if(missing(obj)){
		stop("obj is missing");
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object");
		}
		ncell = nrow(obj);
		if(!isDimReductComplete(obj@smat)){
			stop("dim.reduct is not complete, run 'runDimReduct' first")
		}
	}
		
	down.sample = min(down.sample, ncell);
	idx.ds = sort(sample(seq(ncell), down.sample));
	obj = obj[idx.ds,,drop=FALSE];

	if(max(pca.dims) > length(obj@smat@sdev)){
		stop(paste("pca.dims exceeds PCA dimentions ", length(obj@smat@sdev)));
	}
	
	if((x=length(pca.dims)) > 50L){
		stop("pca.dims must be within 1:50")
	}
	
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	op <- par(mfrow = c(5,5), oma = c(3,3,1,1) + 0.2, mar = c(0,0,1,1) + 0.2);
	PCA.plot <- split(sort(pca.dims), ceiling(seq(pca.dims)/2));
	if((length(x = pca.dims)  %% 2) == 1){
		PCA.plot = PCA.plot[1:(length(PCA.plot) - 1)]
	}
	
	for(x in PCA.plot){
		data.use = data.frame(obj@smat@dmat[,c(x[1],x[2])]);	
		colnames(data.use) = c("dim1", "dim2");
		plot(x=data.use[,1], 
			 y=data.use[,2],
			 cex=point.size, 
			 col=scales::alpha(point.color, point.alpha),
			 mtext(paste(paste("PC", x[1]), x[2], sep=" vs "), side=3),
			 yaxt='n', 
			 xaxt="n",
			 xlab="", 
			 ylab=""
		);
	}

	if(!is.null(pdf.file.name)){
		dev.off()		
	}	
	graphics::par(mfrow=c(1,1));
}

#' Plot gene-body accessibility level
#'
#' @param obj A snap object.
#' @param gene.names Name of genes to plot.
#' @param viz.method Visulization method c("tsne", "umap").
#' @param point.size Point size [0.5].
#' @param point.shape Point shape type [19].
#' @param point.color Point color ["blue"].
#' @param background.point If add points as background [TRUE].
#' @param background.point.color Color of background points ["grey"].
#' @param background.point.alpha Transparency level of background points [0.3].
#' @param background.point.size Size of background points [0.5].
#' @param background.point.shape Shape of background points [19].
#' @param low.value Feature value is standarded to 0-1, value less than low.value will be set to low.value.
#' @param high.value Feature value is standarded to 0-1, value greater than high.value will be set to high.value.
#' @param down.sample Number of cells to plot. Cells will be randomly downsampled for this number for plotting.
#' @param seed.use Random seed [10].
#' @param plot.nrow Number of rows in the plot [3].
#' @param plot.ncol Number of columns in the plot [3].
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @param ... Arguments passed to plot method.
#' 
#' @examples 
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	mat="bmat"
#'	);
#' demo.sp = runNormJaccard(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir()
#'	);
#' demo.sp = runDimReduct(
#'	obj=demo.sp, 
#'	pc.num=10, 
#'	input.mat="jmat"
#'	);
#' demo.sp = runKNN(
#'	obj=demo.sp, 
#'	pca.dims=1:5, 
#'	k=15, 
#'	snn=TRUE, 
#'	save.knn=FALSE
#'	);
#' demo.sp = runCluster(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	louvain.lib="R-igraph"
#'	);
#' demo.sp = runViz(
#'	obj=demo.sp, 
#'	tmp.folder=tempdir(), 
#'	pca.dims=1:5, 
#'	method="Rtsne"
#'	);
#' demo.sp = scaleCountMatrix(
#'	obj=demo.sp,
#'	mat="gmat",
#'	cov=rowSums(demo.sp, mat="bmat"),
#'	method="RPM"
#'	)
#' gene.names = c(
#'	"Prdm14",   "E330040D14Rik", "Gm17971",       
#'	"Gm17970",	"Defb44-ps",     "Gm7357",
#'	"Gm37265",  "Kctd18",		"Gm37143"
#'	)
#' plotGene(
#'	obj=demo.sp, 
#'	gene.names=gene.names, 
#'	viz.method="tsne"
#'	)
#' @importFrom grDevices pdf dev.off
#' @importFrom methods slot 
#' @importFrom scales alpha
#' @importFrom graphics par points plot 
#' @importFrom grDevices pdf dev.off
#' @export
plotGene <- function(obj, gene.names, viz.method, point.size, point.color, point.shape,  background.point, background.point.color, background.point.alpha, background.point.size, background.point.shape, low.value, high.value, down.sample, seed.use, plot.nrow, plot.ncol, pdf.file.name, pdf.height, pdf.width,...){
  UseMethod("plotGene", obj);
}

#' @export
plotGene.default <- function(
	obj, 
	gene.names,
	viz.method=c("tsne", "umap"),
	point.size=0.5,
	point.color="red",
	point.shape=19,
	background.point=TRUE,
	background.point.color="grey",
	background.point.alpha=0.3,
	background.point.size=0.5,
	background.point.shape=19,
	low.value=0.0,
	high.value=1.0,
	down.sample=5000,
	seed.use=10,
	plot.nrow=3,
	plot.ncol=3,
	pdf.file.name=NULL, 
	pdf.height=7, 
	pdf.width=7,
	...
){
	
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
		ncell = nrow(obj);
		data.use = obj@gmat;
		if((x=nrow(data.use)) == 0L){
			stop("gmat is empty, run addGamtToSnap frist")
		}
		if((x = nrow(data.use))!= ncell){
			stop("gmat has different number of rows from cell number, add gmat again")
		}
	}
	
	if(missing(gene.names)){
		stop("gene.names is missing")
	}else{
		if((x=length(gene.names)) > 9L){
			warning("plotGene works best with less than 9 genes in on plot");
		}
		if(any(!(gene.names %in% colnames(data.use)))){
			stop(paste(gene.names[which(!(gene.names %in% colnames(data.use)))], "does not exist in cell-by-gene matrix"))
		};	
	}

	viz.method = match.arg(viz.method);
	viz.use = methods::slot(obj, viz.method);
	
	if((x = nrow(viz.use)) == 0L){
		stop("visulization matrix is empty, run runViz first")
	}
	
	
	if(down.sample < ncell){
		set.seed(seed.use);
		idx = sort(sample(seq(ncell), down.sample));
		data.use = data.use[idx,];
		viz.use = viz.use[idx,];
	}
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	ndim = length(gene.names);
	
	op <- par(mfrow = c(plot.nrow,plot.ncol), oma = c(3,3,1,1) + 0.2, mar = c(0,0,1,1) + 0.2);
	
	for(i in seq(ndim)){
		y = data.use[,which(colnames(data.use) == gene.names[i])];
		y = pmin(1, y/quantile(y[which(y > 0)], 0.99));
		y[y < low.value] = 0;
		y[y > high.value] = 1;
		if(background.point){
			plot(viz.use, 
	   			 main=gene.names[i],
				 col=scales::alpha(background.point.color, background.point.alpha), 
				 cex=background.point.size,
				 pch=background.point.shape,
		   		 yaxt='n', 
		   		 xaxt="n",
		   		 xlab="", 
		   		 ylab="",
				 ...
				 );	
		 	points(viz.use, 
		 		   col=alpha(point.color, y), 
		 		   cex=point.size,
		 		   pch=point.shape
		 		 );	
		}else{
			plot(viz.use, 
	   			 main=gene.names[i],
				 col=alpha(point.color, y), 
				 cex=point.size,
				 pch=point.shape,
		   		 yaxt='n', 
		   		 xaxt="n",
		   		 xlab="", 
		   		 ylab="",
				 ...
				 );	
		}
	}
	
	if(!is.null(pdf.file.name)){
		dev.off()		
	}
	par(mfrow=c(1,1));
}

#' Feature Enrichment Boxplot
#'
#' @param obj A snap object.
#' @param feature Feature enrichment value for each cell.
#' @param outline If 'outline' is not true, the outliers are not drawn (as points whereas S+ uses lines).
#' @param ylab Name of ylab.
#' @param main Main title.
#' @param add.point If 'add.point' is true, scatter points will be added to the top of the boxplot.
#' @param point.size Point size [0.5].
#' @param point.shape Point shape type [19].
#' @param point.alpha Point transparancy level [0.9].
#' @param pdf.file.name pdf file name to save the plot [NULL].
#' @param pdf.width the width of the graphics region in inches [7].
#' @param pdf.height the height of the graphics region in inches [7].
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics stripchart boxplot
#' @importFrom methods slot 
#' @importFrom scales alpha
#' @importFrom methods is
#' @export
boxPlotFeature <- function(obj, feature, outline, ylab, main, add.point,  point.size, point.shape, point.alpha, pdf.file.name, pdf.height, pdf.width){
  UseMethod("boxPlotFeature", obj);
}

#' @export 
boxPlotFeature.default <- function(
	obj, 
	feature, 
	outline=FALSE, 
	ylab=NULL, 
	main=NULL, 
	add.point=TRUE, 
	point.size=0.2, 
	point.shape=19, 
	point.alpha=0.5, 
	pdf.file.name=NULL, 
	pdf.height=7, 
	pdf.width=7
){
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap object")
		}
		ncell = nrow(obj);
		data.use = obj@cluster;
		if((x=length(data.use)) == 0L){
			stop("cluster is empty, run runCluster frist")
		}
		if((x = length(data.use))!= ncell){
			stop("gmat has different number of rows from cell number, add gmat again")
		}
	}
	
	if(length(data.use) != length(feature)){
		stop("feature has different length with cell number");
	}
	
	
	if(!is.null(pdf.file.name)){
		if(file.exists(pdf.file.name)){
			warning("pdf.file already exists");
			file.remove(pdf.file.name);
		}else{
			if(!file.create(pdf.file.name)){
				stop("cannot create pdf.file, not a directory")				
			}
			file.remove(pdf.file.name);
		}	
		grDevices::pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
	}
	
	colPanel = createColorPanel(length(unique(obj@cluster)));	
	
	data.use = data.frame(y=feature, cluster=data.use);
	graphics::boxplot(y ~ cluster, data=data.use, las=3, outline=outline, xlab=NULL, mian=main, ylab=ylab, col = colPanel[data.use$cluster]);
	if(add.point){
		graphics::stripchart(y ~ cluster,  data=data.use, vertical = TRUE, method = "jitter", add = TRUE, pch = point.shape, cex=point.size, col = alpha(colPanel[data.use$cluster], point.alpha))		
	}
	
	if(!is.null(pdf.file.name)){
		grDevices::dev.off()		
	}
	graphics::par(mfrow=c(1,1));
}


######' Feature Enrichment scatterplot
######'
######' @param obj A snap object.
######' @param method Visulziation method c("tsne", "umap").
######' @param feature Feature enrichment value for each cell.
######' @param point.size Point size [0.5].
######' @param point.shape Point shape type [19].
######' @param point.color Point color ["blue"].
######' @param pdf.file.name pdf file name to save the plot [NULL].
######' @param pdf.width the width of the graphics region in inches [7].
######' @param pdf.height the height of the graphics region in inches [7].
######' @param ... Arguments passed to plot function.
######' @importFrom grDevices pdf dev.off
######' @importFrom graphics stripchart boxplot
######' @importFrom methods slot 
######' @importFrom cowplot plot_grid
######' @importFrom scales alpha
######' @export
#####scatterPlotFeature <- function(obj, method, feature, point.size, point.shape, pdf.file.name, pdf.height, pdf.width, ...){
#####  UseMethod("scatterPlotFeature", obj);
#####}
#####
#####scatterPlotFeature.default <- function(
#####	obj,
#####	method=c("tsne", "umap"),
#####	point.size=1,
#####	point.shape=19,
#####	pdf.file.name=NULL,
#####	pdf.height=7,
#####	pdf.width=7,
#####	...
#####){
#####	
#####}
#####


############################################################
###' Project Gene-body Accessibility to T-SNE plot
###'
###' @param obj A snap object.
###' @param gene.names Name of genes to plot.
###' @param viz.method Project method c("tsne", "umap").
###' @param point.size Point size [0.5].
###' @param point.shape Point shape type [19].
###' @param point.color Point color ["blue"].
###' @param point.alpha Point transparancy level [0.9].
###' @param down.sample Number of cells to plot. Cells will be randomly downsampled for this number for plotting.
###' @param pdf.file.name pdf file name to save the plot [NULL].
###' @param pdf.width the width of the graphics region in inches [7].
###' @param pdf.height the height of the graphics region in inches [7].
###' @param ... Arguments passed to ggplot.
###' @importFrom grDevices pdf dev.off
###' @import ggplot2 
###' @importFrom methods slot 
###' @importFrom cowplot plot_grid
###' @export
##plotGene <- function(obj, gene.names, viz.method, point.size, point.color, point.shape, point.alpha, down.sample, pdf.file.name, pdf.height, pdf.width,...){
##  UseMethod("plotGene", obj);
##}
##
###' @export
##plotGene.default <- function(
##	obj, 
##	gene.names,
##	viz.method=c("tsne", "umap"),
##	point.size=0.5,
##	point.color="blue",
##	point.shape=19,
##	point.alpha=0.9,
##	low.value=0.0,
##	high.value=1.0,
##	down.sample=10000,
##	pdf.file.name=NULL, 
##	pdf.height=7, 
##	pdf.width=7,
##	...
##){
##	ncell = nrow(obj);
##	if(missing(obj)){
##		stop("obj is missing");
##	}else{
##		if(!is.snap(obj)){
##			stop("obj is not a snap object")
##		}
##		if((x=nrow(obj@gmat)) == 0L){
##			stop("gmat is empty, run addGamt frist")
##		}
##		data.use = obj@gmat;
##		if((x = nrow(data.use))!= ncell){
##			stop("gmat has different number of rows from cell number")
##		}
##	}
##
##	viz.use = methods::slot(obj, viz.method);
##	if((x = nrow(viz.use))==0L){
##		stop("visulization method does not exist, run runViz first")
##	}
##	
##	if(missing(gene.names)){
##		stop("gene.names is missing")
##	}else{
##		if((x=length(gene.names)) > 9L){
##			warning("plotGene can only plot up to 9 genes in on plot");
##			gene.names = gene.names[1:9];
##		}
##		if(any(!(gene.names %in% colnames(data.use)))){
##			stop(paste(gene.names[which(!(gene.names %in% colnames(data.use)))], "does not exist in cell-by-gene matrix"))
##		};	
##	}
##	
##	if(!is.null(pdf.file.name)){
##		if(file.exists(pdf.file.name)){
##			warning("pdf.file already exists");
##			file.remove(pdf.file.name);
##		}else{
##			if(!file.create(pdf.file.name)){
##				stop("cannot create pdf.file, not a directory")				
##			}
##			file.remove(pdf.file.name);
##		}	
##		pdf(pdf.file.name,width=pdf.width,height=pdf.height); 
##	}
##
##	plt.list = lapply(as.list(seq(length(gene.names))), function(i){
##		idy = which(colnames(data.use)==gene.names[i]);
##		y = data.use[,idy];
##		y[y==0] = NA;
##		y = pmin(1, y / quantile(y, 0.98, na.rm=TRUE));
##		y[y > high.value] = high.value;
##		y[y < low.value] = low.value;
##		data.use = data.frame(viz.use[,1], viz.use[,2], y);
##		colnames(data.use) = c("tsne1", "tsne2", "value");
##		data.use.pos = data.use[which(y > 0),];
##		sp2 <- ggplot() + 
##			geom_point(data=data.use, aes(x=data.use$tsne1, y=data.use$tsne2), color="#D3D3D3", size=point.size) +
##			geom_point(data=data.use.pos, aes(x=data.use.pos$tsne1, y=data.use.pos$tsne2, color= data.use.pos$value), size=point.size, shape=point.shape, alpha=point.alpha) +
##			scale_colour_gradient2() + 
##  		  	theme(
##				panel.grid.major = element_blank(), 
##				panel.grid.minor = element_blank(),
##				panel.background = element_blank(), 
##				legend.position='none', 
##				axis.line = element_line(colour = "black", linetype = "solid")
##				) +			
##			ggtitle(label=gene.names[i]) +
##		    xlab("Dim-1") + ylab("Dim-2")
##		return(sp2)
##	})
##	
##    if (length(x = plt.list) > 9) {
##      nCol <- 4
##    } else {
##      nCol <- min(length(x = plt.list), 3)
##    }
##	
##	plots.combined <- plot_grid(plotlist = plt.list, ncol = nCol);
##	if(!is.null(pdf.file.name)){
##		dev.off()		
##	}else{
##		print(plots.combined)
##	}	
##}
