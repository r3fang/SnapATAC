#' @export
fastMatMult <- function(A, B){
	if(length(find("eigenMapMatMult")) == 0){
		sourceCpp(code='
			// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
			#include <RcppArmadillo.h>
			#include <RcppEigen.h>
	
			typedef Eigen::MappedSparseMatrix< double > mappedSparseMatrix ;
			typedef Eigen::Map< Eigen::VectorXd > mappedVector ;
	
			// [[Rcpp::export]]
			SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
			    Eigen::MatrixXd C = A * B;
			    return Rcpp::wrap(C);
			}'
		)		
	}
	return(eigenMapMatMult(A, B));
} 

