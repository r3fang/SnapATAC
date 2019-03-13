####### Add Halo to Text
####### @param plt a ggplot object
####### @param x,y numeric vectors of coordinates where the text 'labels' should be written.  If the length of 'x' and 'y' differs, the shorter one is recycled.
####### @param labels a character vector or expression specifying the _text_ to be written.  An attempt is made to coerce other language objects (names and calls) to expressions, and vectors and other classed objects to character vectors by 'as.character'.  If 'labels' is longer than 'x' and 'y', the coordinates are recycled to the length of 'labels'.
####### @param text.color text color to use ["black"]
####### @param text.size text size to use [1]
####### @param text.halo To add halo around text [TRUE].
####### @param text.halo.color Halo color ["white"].
####### @param text.halo.width Halo with [0.2].
######textHalo_ggplot2 <- function(
######	plt,
######	x, y = NULL, 
######	labels = NULL, 
######	text.color='black', 
######	text.size=1, 
######	text.halo=TRUE, 
######	text.halo.color='white', 
######	text.halo.width=0.1
######){
######	theta <- seq(0, 2*pi, length.out=50);
######	xo <- text.halo.width * diff(range(x))/200
######	yo <- text.halo.width * diff(range(y))/200
######	
######	if(text.halo){
######		for(i in theta) {
######		    plt <- plt + annotate("text", x=x+ cos(i)*xo,y=y+ sin(i)*yo, label=labels, size=text.size, colour=text.halo.color)
######		}
######	}
######	plt <- plt + annotate("text", x=x,y=y,label=labels, size=text.size, colour=text.color)
######	return(plt);
######}

#' @importFrom grDevices xy.coords 
#' @importFrom graphics strwidth strheight 
textHalo <- function(
	x, y=NULL, 
	labels, 
	col='white', 
	bg='black', 
	r=0.1,
	... 
){

	theta= seq(0, 2*pi, length.out=50);
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}


# Create a color panel
# @param num.color Number of diffenrent colors to return
createColorPanel <- function(num.color){
	colPanel = c(
		"#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
		"#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744", 
		"#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
	    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
	    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
	    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
	    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
	    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
	    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
	    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
	    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
	    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
	    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
	   )
	if(num.color > length(colPanel)){
		colPanel = c(colPanel, colVector(num.color - length(colPanel)));
	}else{
		colPanel = colPanel[1:num.color];
	}
	return(colPanel)
}


# Create a color panel
# @param num.color Number of diffenrent colors to return
# @param type Type of colors to generate
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
colVector <- function(num.color=60, type=c("qual", "div", "seq")){
	type <- match.arg(type);
	set.seed(10)
	qual_col_pals = brewer.pal.info[brewer.pal.info$category == type,];
	set.seed(10)
	col_vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))));
	return(col_vector)
} 





