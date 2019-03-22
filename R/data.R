#' Single Nucleus ATAC-seq Object
#'
#' Data from single nucleus ATAC-seq experiment on secondary motor cortex
#' in adult mouse brain. The data is down-sampled to 100 cells with 1,000 bins,
#' 1,000 peaks and 1,000 genes from the original dataset. 
#'
#' @docType data
#'
#' @usage data(demo)
#'
#' @format An object of class \code{"snap"} with 100 cells:
#' \describe{
#'   \item{metData}{meta data of each cell}
#'   \item{bmat}{cell x bin matrix}
#'   \item{pmat}{cell x peak matrix}
#'   \item{gmat}{cell x gene matrix}
#'   \item{jmat}{jaccard index matrix object}
#'   \item{graph}{knn graph}
#'   \item{cluster}{cluster result}
#'   \item{tsne}{tsne cooridnates}
#'   \item{umap}{umap cooridnates}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(demo)
"demo"
