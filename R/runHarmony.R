#' Perform batch effect using harmony
#'
#' Integration of single cell ATAC-seq data for batch correction using Harmony algorithm
#'
#' @param obj A snap object
#' @param eigs.dim A vector of the components to use in construction of the KNN graph.
#' @param meta_data Either (1) Dataframe with variables to integrate or (2) vector with labels.
#' @param vars_use If meta_data is dataframe, this defined which variable(s) to remove (character vector).
#' @param ... Paramters passed to Harmony
#' 
#' @import Matrix
#' @export
runHarmony <- function(obj, eigs.dim, meta_data, vars_use, ...){
	UseMethod("runHarmony", obj);
}

#' @export
runHarmony.default <- function(obj, eigs.dim, meta_data, vars_use=NULL, ...){	
	if (!requireNamespace("harmony", quietly = TRUE)) {
	      stop("Please install harmony first - learn more at https://github.com/immunogenomics/harmony")
	}else{
		require(harmony);		
	}

	cat("Epoch: checking input parameters\n", file = stderr())
	if(missing(obj)){
		stop("obj is missing")
	}else{
		if(!is(obj, "snap")){
			stop("obj is not a snap obj")
		}
	}
	
	if(!(isDimReductComplete(obj@smat))){
		stop("dimentionality reduction is not complete, run 'runDiffusionMaps' first")
	}

	ncell = nrow(obj);
	nvar = dimReductDim(obj@smat);
	
	if(missing(eigs.dim)){
		stop("eigs.dim is missing")
	}else{
		if(is.null(eigs.dim)){
			eigs.dim=1:nvar;	
		}else{
			if(any(eigs.dim > nvar) ){
				stop("'eigs.dim' exceeds component number");
			}		
		}
	}
		
	data.use = weightDimReduct(obj@smat, eigs.dim, FALSE);

	my_harmony_embeddings <- harmony::HarmonyMatrix(data.use, meta_data=meta_data, vars_use=vars_use, do_pca=FALSE);
	obj@smat@dmat <- my_harmony_embeddings;
	obj@smat@sdev <- obj@smat@sdev[eigs.dim];
	return(obj);
}

