###############################################################################!
#' .checkmateSingleMetaIntegratorObject
#'
#' @description
#' A function with multiple checkmate assertions that verify if an input
#' object is a MetaIntegrator object that contains only a single dataset.
#'
#' @details
#' MetaIntegrator objects aren't a specific class, so can't perform a single
#' class check. Instead they represent a list that follows a standard pattern.
#' This function checks several of these expected patterns.
#'
#' @keywords geo_object_helpers
#'
#'
#' @param geo_object List of format created by MetaIntegrator::getGEOData
#'
#' @return NULL or fails
#'
#' @examples
#' \dontrun{
#'   .checkmateSingleMetaIntegratorObject(GSE6629_GEO) # Returns NULL
#'   .checkmateSingleMetaIntegratorObject(GSE6629_expressionRaw) # Returns error
#' }

.checkmateSingleMetaIntegratorObject <- function(geo_object){
  checkmate::assertNames(names(geo_object), must.include = "originalData",
                         .var.name = "Must be a named list with an 'originalData' slot")
  checkmate::assertList(geo_object$originalData, len = 1,
                        .var.name = "Must only have one GSE slot under 'originalData'")
  checkmate::assertTRUE(
    MetaIntegrator::checkDataObject(
      geo_object$originalData[[1]],
      objectType = "Dataset"),
    .var.name = "GSE must be a MetaIntegrator dataset")
  # checkmate::assertNames(names(geo_object$originalData[[1]]),
  #                        must.include = "expr",
  #                        .var.name = "GSE must include an 'expr' slot")
}

###############################################################################!
#' addRawExprMatrix
#'
#' @description
#' Add an exprRaw slot to a geo_object containing a user specified
#' raw data expression matrix
#'
#' @keywords geo_object_helpers
#'
#' @param geo_object List of format created by MetaIntegrator::getGEOData
#' @param expr Matrix containing raw expression data
#'
#' @return List of format created by getGEO with additional $exprRaw slot
#' @export
#'
#' @examples
#' data(GSE6629_GEO)
#' data(GSE6629_expressionRaw)
#' gse <- addRawExprMatrix(GSE6629_GEO, GSE6629_expressionRaw)
#' "exprRaw" %in% names(gse$originalData[[1]]) # TRUE
#' gse$originalData[[1]]$exprRaw == GSE6629_expressionRaw # TRUE
addRawExprMatrix <- function(geo_object, expr, mode = "microarray"){

  # Input checks
  .checkmateSingleMetaIntegratorObject(geo_object)

  # Other checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assert_choice(mode, c("microarray", "rnaseq"), add = all_checks)
  if (mode == "microarray"){
    checkmate::assert_matrix(expr,
                             ncols = ncol(geo_object$originalData[[1]]$expr),
                             min.rows = nrow(geo_object$originalData[[1]]$expr),
                             add = all_checks)
  }

  if (mode == "rnaseq"){
    checkmate::assert_matrix(expr, add = all_checks)
    # TODO Add better check. Instead of above, since RNAseq won't have an exp matrix, make it from pheno instead?
  }

  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  # Update geo_object
  geo_object$originalData[[1]]$exprRaw <- expr

  return(geo_object)
}

###############################################################################!
#' addReprocessedExprMatrix
#'
#' @description
#' Update the $expr slot in a geo_object with a new expression matrix
#'
#' @details
#' getGEO returns an object with an expression matrix processed according to
#' the original study investigators. However, users may want to perform their
#' own processing of raw data (e.g. standardization for microarray CEL files).
#' This function updates a geo_object $expr matrix with a user provided matrix
#' and reassigns the original, author-provided matrix to $exprGEO. Adds
#' additional explanatory comments to $expr_comments and $exprGEO_comments.
#' $ Uses `checkmate` to assert certain conditions.
#'
#' @keywords geo_object_helpers
#'
#' @param geo_object List of format created by getGEO
#' @param expr Matrix. Matching dimensions of $expr in geo_object and rownames
#' all found within names of $keys (i.e. probes) or geo_object
#'
#' @return List of format created by getGEO with additional $exprGEO and
#' $expreGEO_comment slots
#' @export
#'
#' @examples
#' gse <- addReprocessedExprMatrix(GSE6629_GEO, GSE6629_expressionNormalized)
#' "exprGEO" %in% names(gse$originalData[[1]]) # TRUE
#' gse$originalData[[1]]$exprGEO == GSE6629_GEO$originalData[[1]]$expr # TRUE
addReprocessedExprMatrix <- function(geo_object, expr, mode = "microarray") {

  # Check geo_object
  .checkmateSingleMetaIntegratorObject(geo_object)

  # Additional checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assert_choice(mode, c("microarray", "rnaseq"), add = all_checks)

  ## Check replacement matrix properties
  ### Only if microarray, since RNAseq GEO objects won't have an original matrix
  if (mode == "microarray"){
    original_expr <- geo_object$originalData[[1]]$expr
    original_comment <- geo_object$originalData[[1]]$expr_comment
    probe_names <- names(geo_object$originalData[[1]]$keys)
    dim_check <- all(dim(original_expr) == dim(expr))
    probe_check <- all(rownames(expr) %in% probe_names)

    checkmate::assert_true(dim_check, add = all_checks,
                           .var.name = "Replacement expr dimensions must match original expr dimensions")
    checkmate::assert_true(probe_check, add = all_checks,
                           .var.name = "Replacement expr probes (i.e. rownames) must all be found in original $keys")
  }

  if (mode == "rnaseq"){
    original_expr <- NULL
    original_comment <- NULL
  }

  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  # Update geo_object
  geo_object$originalData[[1]]$exprGEO <- original_expr
  geo_object$originalData[[1]]$exprGEO_comment <- original_comment
  geo_object$originalData[[1]]$expr <- expr
  geo_object$originalData[[1]]$expr_comment <-
    "$expr replaced with reprocessed expression matrix from raw data.
  Original GEO expression matrix moved to $exprGEO and $exprGEO_comment"

  return(geo_object)
}

###############################################################################!

#' fixMetaIntegratorForRNAseq
#'
#' @description
#' Add minimal necessary empty-data to a MetaIntegrator object created from an
#' RNAseq dataset in order to pass `MetaIntegrator::checkDataObject`
#'
#' @details
#' `MetaIntegrator::checkDataObject` is a convenient function to check that an
#' object is a valid `MetaIntegrator` object. However, it was designed with
#' microarray data in mind. When a microarray dataset is downloaded from GEO
#' using `MetaIntegrator`, an expression matrix is always included (`expr` slot)
#' as well as the probe-gene mappings (`keys` slot). However, for RNAseq data,
#' GEO does not include a standard expression matrix, this must be downloaded
#' specifically from supplemental files. In addition, there are no probe
#' mappings, since RNAseq doesn't use probes.
#'
#' This function creates near-empty data in the expected slots to all
#' a MetaIntegrator object on RNA seq data to pass the
#' `MetaIntegrator::checkDataObject` check
#'
#' @keywords geo_object_helpers
#'
#' @param geo_object
#'
#' @return MetaIntegrator object
#' @export
#'
#' @examples
#' gse_rna_original <- GSE158395_GEO
#' gse_rna_fixed <- fixMetaIntegratorForRNAseq(GSE158395_GEO)
#' # Returns FALSE with multiple warnings
#' MetaIntegrator::checkDataObject(gse_rna_original$originalData[[1]], "Dataset")
#' # Returns TRUE
#' MetaIntegrator::checkDataObject(gse_rna_fixed$originalData[[1]], "Dataset")
fixMetaIntegratorForRNAseq <- function(geo_object){
  geo_object$originalData[[1]]$keys <-  c("rnaseq" = "RNASEQ")
  geo_object$originalData[[1]]$expr <- matrix(c(0))
  rownames(geo_object$originalData[[1]]$expr) <- "rnaseq"
  return(geo_object)
}
