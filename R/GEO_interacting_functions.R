

#' isValidGSE
#'
#' Check whether a GEO accession is valid
#'
#' @param accession String of GEO accession to test
#'
#' @return Logical
#' @export
#'
#' @examples
#' \donttest{
#' isValidGSE("GSE18606") # Return TRUE
#' isValidGSE("GSE18606a") # Return FALSE
#' }
isValidGSE <- function(accession) {
  # Accession must start with GEO
  if (!grepl("^GSE", accession)){return(FALSE)}
  # Accession must have an entry
  search_result <- rentrez::entrez_search(db = "gds", term = accession)
  return(length(search_result$ids) > 0)
}
