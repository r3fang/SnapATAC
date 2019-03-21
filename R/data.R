#' Single Nucleus ATAC-seq for MOs 2K
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
#' data(mos)
"mos"
