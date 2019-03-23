#' Single Nucleus ATAC-seq Object
#'
#' Data from single nucleus ATAC-seq experiment on secondary motor cortex
#' in adult mouse brain. The data is down-sampled to 171 cells with 27,348 100kb bins,
#' 999 peaks and 1,000 genes. 
#'
#' @docType data
#'
#' @usage data(demo.sp)
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
#' data(demo.sp)
"demo.sp"
