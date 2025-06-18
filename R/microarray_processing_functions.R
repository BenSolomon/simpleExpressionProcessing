###############################################################################!
#' processAffy
#'
#' @description
#' Processing pipeline for raw affymetrics array data. Utilizes affy package
#'
#' @param cel_dir String. Directory containing CEL files to process
#'
#' @return List of matrices. Two expression matrices from processed affy
#' object, one of raw data, one of normalized data
#' @export
#'
#' @examples
#' \donttest{
#' geo_accession <- "GSE6629"
#' temp_dir <- getwd()
#' gse_dir <- sprintf("%s/%s", temp_dir, geo_accession)
#'
#' getSuppFiles(GEO = geo_accession, #TODO Suppress curl download message
#'              makeDirectory = TRUE,
#'              baseDir = temp_dir,
#'              fetch_files = TRUE,
#'              filter_regex = NULL,
#'              download_method = "curl",
#'              unpack = TRUE)
#'
#' affy_output_list <- processAffy(gse_dir)
#' names(affy_output_list) # c("raw_expression", "normalized_expression")
#' dim(affy_output_list$raw_expression) # c(1354896, 2)
#' dim(affy_output_list$normalized_expression) # c(54675, 2)
#' }
processAffy <- function(cel_dir){
  checkmate::assertDirectoryExists(cel_dir)
  study <- basename(cel_dir)
  cel_files <- affy::list.celfiles(cel_dir, full.names = T)
  affy_array <- affy::ReadAffy(filenames = cel_files)

  Biobase::sampleNames(affy_array) <- stringr::str_extract(
    Biobase::sampleNames(affy_array),
    "^GSM[0-9]+"
  )

  raw_expression <- affy::exprs(affy_array)
  affy_array <- affy::rma(affy_array)
  rma_expression <- affy::exprs(affy_array)

  output_list <- list(
    "raw_expression" = raw_expression,
    "normalized_expression" = rma_expression
  )

  return(output_list)
}
