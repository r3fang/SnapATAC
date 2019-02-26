#' Single Nucleus ATAC-seq for MOs
#'
#' Data from single nucleus ATAC-seq experiment on secondary motor cortex
#' in adult mouse brain. The data is down-sampled to 2,000 cells from the 
#' original dataset. 
#'
#' @docType data
#'
#' @usage data(mos)
#'
#' @format An object of class \code{"snap"} with 2000 cells:
#' \describe{
#'   \item{metData}{meta data of each cell}
#'   \item{bmat}{cell x bin matrix}
#'   \item{pmat}{cell x peak matrix}
#'   \item{gmat}{cell x gene matrix}
#'   \item{jmat}{jaccard index matrix}
#'   \item{nmat}{normlaized jaccard index matrix}
#'   \item{cluster}{cluster result}
#'   \item{tsne}{tsne cooridnates}
#'   \item{umap}{umap cooridnates}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(mos)
#' mos
#' number of barcodes: 2000
#' number of bins: 476595
#' number of peaks: 316257
#' number of genes: 53278
#' ==========================
#' meta data            (metData) :  TRUE
#' cellxbin matrix      (bmat)    :  TRUE
#' cellxpeak matrix     (pmat)    :  TRUE
#' cellxgene matrix     (gmat)    :  TRUE
#' jaccard matrix       (jmat)    :  FALSE
#' normalization        (nmat)    :  FALSE
#' PCA:                 (smat)    :  FALSE
#' cluster:             (cluster) :  FALSE
#' t-sne:               (tsne)    :  FALSE
#' umap:                (umap)    :  FALSE
"mos"
