#' geoAllBoxplots
#'
#' @description
#' Starting with a MetaIntegrator GEO object, plots a boxplot of expression
#' values for all samples in an expression #' matrix, broken up by which
#' expression matrix, and which dataset.
#'
#' @param geo_object List of format created by MetaIntegrator::getGEOData
#'
#' @return Base R plot
#' @export
#'
#' @examples
#' \donttest{
#'   # Will plot a single boxplot for the original expression matrix
#'   geoAllBoxplots(GSE6629_GEO)
#' }
geoAllBoxplots <- function(geo_object){
  # TODO Create a variation of .checkmateSingleMetaIntegratorObject, in order
  # to validate that geo_object is a valid MetaIntegrator object

  # For all datasets in geo_object
  geo_accessions <- names(geo_object$originalData)
  for (g in geo_accessions){
    # Get all slot names in dataset
    slots <- names(geo_object$originalData[[g]])
    # Evaluate if each slot is an "expr" slot
    expr_slot <- grepl("^expr", slots)
    # Evaluate if each slot is a matrix
    matrix_slot <- unname(sapply(geo_object$originalData[[g]], is.matrix))
    # Combine to determine which slots are expr matrices
    # Double check necessary to exclude `expr comment` slot, if present
    expr_matrices <- slots[expr_slot & matrix_slot]
    # For all expression matrices in dataset
    for (e in expr_matrices){
      title <- sprintf("%s - %s", g, e)
      boxplot(geo_object$originalData[[g]][[e]], outline = F, main = title)
    }
  }
}

