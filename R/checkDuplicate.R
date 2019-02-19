#' Check & Remove Duplicates
#'
#' This function tables a data matrix as input and check weather duplicates exist. If so, add white noise to the duplicates
#'
#' @param x the numeric matrix that provides data.
#' @param sd standard deviation for the noise.
#'
#' @export

checkDuplicate <- function(x, ...) {
  UseMethod("checkDuplicate");
}

checkDuplicate.default <- function(x, sd=1e-10){
	dup_id <- which(duplicated(x))
	if (length(dup_id) > 0) {
	       x[dup_id, ] <- x[dup_id, ] + rnorm(length(x) * ncol(x), sd)
	}
	return(x)
}
