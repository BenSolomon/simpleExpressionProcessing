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
#' # Download an Affymetrix dataset
#' geo_accession <- "GSE6629"
#' temp_dir <- getwd()
#' gse_dir <- sprintf("%s/%s", temp_dir, geo_accession)
#' getSuppFiles(GEO = geo_accession,
#'              makeDirectory = TRUE,
#'              baseDir = temp_dir,
#'              fetch_files = TRUE,
#'              filter_regex = NULL,
#'              download_method = "curl",
#'              unpack = TRUE)
#'
#' # Process dataset
#' affy_output_list <- processAffy(gse_dir)
#'
#' # Test output
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

###############################################################################!
#' processAgilent
#'
#' @description
#' Processing pipeline for raw agilent array data. Utilizes limma package
#'
#' @details
#' Unlike raw affy files, raw agilent files do not contain info about probe
#' mapping, so these have to be obtained manually. Here we choose to utilize
#' the probe mappings that are obtained from GEOquery since these should be
#' kept up to date.
#'
#' This workflow also includes a filter that removes probes if they are
#' not expressed across a set proportion of microarrays. The affy processing
#' pipeline does not have this, so you may chose to set this parameter to 0
#' so that no filtering is performed.
#'
#' There is another filter step that removes probes that match to no gene
#' symbol. However, this is currently commented out since this is not done
#' in the affy workflow and significantly changes the appearance of diagnostic
#' normalization box-plots.
#'
#' @param txt_dir String where raw agilent txt files are located
#' @param p_array_threshold Dbl between 0-1. Represents the proportion of
#' total arrays a probe needs to be expressed within in order to pass
#' filtering. If NULL, will not filter
#' @param filter_nonGene_probes Logical of whether to filter out probes that
#' have no corresponding gene. Default is FALSE since affy pipeline does not
#' do this.
#'
#' @return List of matrices. Two expression matrices from processed affy
#' object, one of raw data, one of normalized data
#' @export
#'
#' @examples
#' \dontrun{
#' # Download an Agilent dataset
#' geo_accession <- "GSE16999"
#' temp_dir <- getwd()
#' gse_dir <- sprintf("%s/%s", temp_dir, geo_accession)
#' getSuppFiles(GEO = geo_accession,
#'              makeDirectory = TRUE,
#'              baseDir = temp_dir,
#'              fetch_files = TRUE,
#'              filter_regex = NULL,
#'              download_method = "curl",
#'              unpack = TRUE)
#'
#' # Process dataset
#' agilent_output_list <- processAgilent(gse_dir, p_array_threshold = 0.5)
#'
#' # Test output
#' names(agilent_output_list) # c("raw_expression", "normalized_expression")
#' dim(agilent_output_list$raw_expression) # c(45015, 2)
#' dim(agilent_output_list$normalized_expression) # c(33164, 2)
#' }
processAgilent <- function(txt_dir,
                           p_array_threshold = 0.5,
                           filter_nonGene_probes = FALSE){

  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertDirectory(txt_dir, add = all_checks)
  study <- basename(txt_dir)
  checkmate::assertTRUE(isValidGSE(study), add = all_checks)
  checkmate::assertDouble(p_array_threshold, lower = 0, upper = 1,
                          null.ok = TRUE, add = all_checks)
  checkmate::assertLogical(filter_nonGene_probes, add = all_checks)
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  # Load files
  txt_files <- list.files(txt_dir, pattern = "GSM.*\\.txt$|\\.txt\\.gz$")
  agilent_array <- limma::read.maimages(txt_files, source = "agilent",
                                        path = txt_dir, green.only = T,
                                        other.columns = "gIsWellAboveBG")
  ## Rename columns
  colnames(agilent_array) <- stringr::str_extract(
    colnames(agilent_array), pattern = "^GSM[0-9]+")

  # Get probe mapping
  ## Does a geo object already exist, if not, get it
  geo_path <- sprintf("%s/../%s.RDS", txt_dir, study)
  if (file.exists(geo_path)){
    gse <- readRDS(geo_path)
  } else {
    gse <- MetaIntegrator::getGEOData(study)
  }
  keys <- gse$originalData[[study]]$keys
  keys <- tibble::enframe(keys, "ProbeName", "symbol")
  keys$inGEO <- TRUE
  ## Join probe mapping from geo obj to genes slot
  agilent_array$genes <- agilent_array$genes %>%
    dplyr::left_join(keys)

  # Store raw expression before normalization
  raw_expression <- agilent_array$E
  rownames(raw_expression) <- agilent_array$genes$ProbeName

  # Normalize data
  agilent_array <- limma::backgroundCorrect(agilent_array, method = "normexp")
  agilent_array <- limma::normalizeBetweenArrays(agilent_array, method = "quantile")

  # Filter data
  ## Exclude control probes and probes in GEO keys
  Control <- agilent_array$genes$ControlType==1L
  InGEO <- agilent_array$genes$inGEO
  agilent_array <- agilent_array[!Control & InGEO, ]
  ## Filter based on probe expression in p_array_threshold proportion of arrays
  ### n_arrays needed = total arrays * proportion of arrays specified, rounded
  ### n_arrays set to a minimum of at least 2
  if (!is.null(p_array_threshold)){
    n_array_threshold <- floor(ncol(agilent_array)*p_array_threshold)
    n_array_threshold <- max(2, n_array_threshold)
    IsExpr <- rowSums(agilent_array$other$gIsWellAboveBG > 0) >= n_array_threshold
    agilent_array <- agilent_array[IsExpr, ]
  }
  ## Filter if probe doesn't have a corresponding gene symbol
  ### Default to FALSE since not part of `affy` pipeline
  if (filter_nonGene_probes){
    NoSymbol <- is.na(agilent_array$genes$symbol)
    agilent_array <- agilent_array[!NoSymbol, ]
  }

  # Final output
  normalized_expression <- agilent_array$E
  rownames(normalized_expression) <- agilent_array$genes$ProbeName
  output_list <- list(
    "raw_expression" = raw_expression,
    "normalized_expression" = normalized_expression
  )
  return(output_list)
}
