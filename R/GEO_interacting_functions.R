###############################################################################!
#' isValidGSE
#'
#' @description
#' Check whether a GEO accession is valid
#'
#' @keywords geo_interacting_functions
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
  checkmate::assertString(accession)
  # Accession must start with GEO
  if (!grepl("^GSE", accession)){return(FALSE)}
  # Accession must have an entry
  search_result <- rentrez::entrez_search(db = "gds", term = accession)
  return(length(search_result$ids) > 0)
}

###############################################################################!
#' isValidGPL
#'
#' @description
#' Check whether or not a GPL platform accession is valid
#'
#' @keywords geo_interacting_functions
#'
#' @param accession String of GPL accession to test
#'
#' @return Logical
#' @export
#'
#' @examples
#' \donttest{
#' isValidGPL("GPL6480") # Return TRUE
#' isValidGPL("GPL648000000") # Return FALSE
#' isValidGPL("GSE18606") # Return FALSE
#' }
isValidGPL <- function(accession) {
  checkmate::assertString(accession)
  # Accession must start with GPL
  if (!grepl("^GPL", accession)){return(FALSE)}
  # Accession must have an entry
  search_result <- rentrez::entrez_search(db = "gds", term = accession)
  return(length(search_result$ids) > 0)
}

###############################################################################!
#' getGPLfromGSE
#'
#' @description
#' Quickly obtain the GPL accession for the platform used in a given
#' GSE code experiment
#'
#' @details
#' GPL accessions can also be obtained from a GEOquery object. However, that
#' requires either downloading the whole dataset or loading a saved dataset.
#' This allows querying of the GPL accession without downloading any data.
#'
#' @keywords geo_interacting_functions
#'
#' @param geo_accession String of GEO accession
#'
#' @return string of GPL accession
#' @export
#'
#' @examples
#' \donttest{
#' getGPLfromGSE("GSE18606") # Return GPL6480
#' }

getGPLfromGSE <- function(geo_accession){
  # Check input
  checkmate::assertTRUE(isValidGSE(geo_accession))
  # Search for GPL within geo_accession
  search_result <- rentrez::entrez_search(db = "gds", term = geo_accession)
  search_result <- rentrez::entrez_summary(db = "gds", search_result$ids[1])
  search_result <- search_result$gpl
  search_result <- sprintf("GPL%s", search_result)
  return(search_result)
}

###############################################################################!
#' parseGPLmetadata
#'
#' @description
#' Parses several metadata fields from the geo webpage entry for a
#' given GPL accession. Used for the purspose of performing a microarray
#' platform check prior to running microarray normalization since
#' different platforms require different methods.
#'
#' @details
#' Somewhat slow for the task. Could also be performed by getGEO in geoQuery
#' package, but this is similarly slow and also involves download of the
#' entire GPL data package in GEO, which is typically many MBs in size. #'
#'
#' @keywords geo_interacting_functions
#'
#' @importFrom magrittr %>%
#'
#' @param gpl_accession String of GPL accession. Validity will be checked by
#' isValidGPL
#'
#' @return data.frame with meta data parameters and values
#' @export
#'
#' @examples
#' \dontrun{
#' parseGPLmetadata("GPL6480") # Return data.frame
#' parseGPLmetadata("GPL648000000") # Error
#' parseGPLmetadata("GSE18606") # Error
#' }
parseGPLmetadata <- function(gpl_accession){
  # Input checks
  checkmate::assertTRUE(isValidGPL(gpl_accession))

  # Construct URL
  url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s", gpl_accession)

  # Parse HTML
  html_doc <- rvest::read_html(url)
  param_value_rows <- html_doc %>%
    rvest::html_nodes("tr") %>%
    .[sapply(., function(row) length(rvest::html_nodes(row, "td")) == 2)]
  parameters <- param_value_rows %>%
    rvest::html_nodes("td:nth-child(1)") %>%
    rvest::html_text(trim = TRUE)
  values <- param_value_rows %>%
    rvest::html_nodes("td:nth-child(2)") %>%
    rvest::html_text(trim = TRUE)

  # Assemble data.frame
  result_df = data.frame(parameter = parameters,
                         value = values,
                         stringsAsFactors = F)
  allowed_parameters <- c("Title", "Organism", "Manufacturer", "Description",
                          "Submission date", "Last update date", "Organization",
                          "Technology type")
  result_df <- result_df[result_df$parameter %in% allowed_parameters, ]
  if (nrow(result_df) == 0){stop("Invalid gpl_accession")}
  return(result_df)
}

###############################################################################!
#' platformCheck
#'
#' @description
#' Attempts to check verify if the platform of a GPL accession matches
#' a query platform
#'
#' @details
#' Intended for uses in input checking for microarray processing. Want to be
#' able to detect
#'
#' @keywords geo_interacting_functions
#'
#' @param gpl_accession String of GPL accession. Validity will be checked by
#' @param platform_query String representing platform of interest
#' @param quiet Logical. Whether to how function is searching for platform
#'
#' @return Logical
#' @export
#'
#' @examples
#'  \donttest{
#' platformCheck("GPL6480", "agilent") # Returns TRUE
#' platformCheck("GPL6480", "affy") # Returns FALSE
#' }
platformCheck <- function(gpl_accession, platform_query, quiet = F){
  df <- parseGPLmetadata(gpl_accession)
  platform_title <- df$value[df$parameter == "Title"]
  if (!quiet){
    message(sprintf("Searching for '%s' within '%s'",
                    platform_query,
                    platform_title))}
  search_result <- grepl(platform_query, platform_title, ignore.case = T)
  return(search_result)
}

###############################################################################!
#' platformGuess

#' @description
#' Attempts to check parese and return the platform of a GPL accession matches
#' a query platform
#'
#' @details
#' Intended for uses in input checking for microarray processing. Want to be
#' able to detect
#'
#' @keywords geo_interacting_functions
#'
#' @param gpl_accession String of GPL accession. Validity will be checked by
#' @param quiet Logical. Whether to how function is searching for platform
#'
#' @return Logical
#' @export
#'
#' @examples
#' \donttest{
#' platformGuess("GPL6480") # Return "agilent"
#' }
platformGuess <- function(gpl_accession, quiet = F){
  patterns <- c("affy", "agilent")
  df <- parseGPLmetadata(gpl_accession)
  platform_title <- df$value[df$parameter == "Title"]
  output_check <- sapply(patterns, function(p) grepl(p, platform_title, ignore.case = T))
  output <- patterns[output_check]
  output <- ifelse(length(output) != 1,
                   NA, output)
  return(output)
}

###############################################################################!
#' getSuppfiles
#'
#' @description
#' Wrapper around GEOquery::getGEOSuppFiles intended to download raw
#' microarray CEL files
#'
#' @keywords geo_interacting_functions
#'
#' @param GEO String. GEO accession
#' @param makeDirectory Logical. Whether to download files to a subfolder for the
#' GEO accession
#' @param baseDir String. Directory to download files or create subdirectory
#' @param fetch_files Logical. Whether to actually download files or return
#' names instead
#' @param filter_regex String. Pattern to filter potential download files for
#' @param download_method String. Takes 'auto', 'wget', libcurl', 'curl'.
#' Specifies the download protocol for getGEOSuppFiles. On my setup, default
#' auto (wget) had connection issues that curl solved
#' @param unpack Logical. Whether to unpack the supplemental tar file that
#' contains the CEL files
#'
#' @return NULL
#' @export
#'
#' @examples
#'  \donttest{
#' geo_accession <- "GSE6629"
#' baseDir <- getwd()
#' getSuppFiles(GEO = geo_accession,
#'              makeDirectory = TRUE,
#'              baseDir = baseDir,
#'              fetch_files = TRUE,
#'              filter_regex = NULL,
#'              download_method = "curl",
#'              unpack = TRUE)
#' output_tar <- sprintf("%s/%s/GSE6629_RAW.tar", baseDir, geo_accession)
#' untar_files <- c("GSM153779.CEL.gz", "GSM153780.CEL.gz")
#' output_untar <- sprintf("%s/%s/%s", baseDir, geo_accession, untar_files)
#' file.exists(output_tar) # TRUE
#' all(file.exists(untar_files)) # TRUE
#' }
getSuppFiles <- function(GEO,
                         makeDirectory = TRUE,
                         baseDir = getwd(),
                         fetch_files = TRUE,
                         filter_regex = NULL,
                         download_method = "curl",
                         unpack = TRUE) {
  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertTRUE(isValidGSE(GEO), add = all_checks)
  checkmate::assertDirectoryExists(baseDir, add = all_checks)
  checkmate::assertLogical(makeDirectory, add = all_checks)
  checkmate::assertLogical(fetch_files, add = all_checks)
  checkmate::assertLogical(unpack, add = all_checks)
  checkmate::assertString(filter_regex, null.ok = TRUE, add = all_checks)
  checkmate::assertChoice(
    download_method,
    choices = c("curl", "auto", "wget", "libcurl"),
    add = all_checks)
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  # Set download method
  options('download.file.method.GEOquery' = download_method) # Auto, wget, and libcurl are slow for me
  # Download supplemental files
  GEOquery::getGEOSuppFiles(
    GEO = GEO,
    makeDirectory = makeDirectory,
    baseDir = baseDir,
    fetch_files = fetch_files,
    filter_regex = filter_regex
  )
  # Unpack tar
  if (unpack) {
    tardir <- sprintf("%s/%s", baseDir, GEO)
    downloaded_files <- list.files(tardir, full.names = T)
    tarfile <- downloaded_files[grepl("tar$", downloaded_files)]
    lapply(tarfile, function(f) untar(tarfile = f, exdir = tardir))
  }
}

###############################################################################!
#' getGEOData_retryWrapper
#'
#' @description
#' Wrapper around getGEOData that will repeat getGEOData for a given accession
#' until data is successfully downloaded, up to n-number of total attempts
#'
#' @details
#' getGEOData can fail, returning 'Warning: Unable to correctly download
#' geo_accession . Dataset will be excluded from this object.' However, this can
#' be due to a connection error, not because no data exists
#'
#' @keywords geo_interacting_functions
#'
#' @param geo_accession GEO accession ID
#' @param n Max number of times getGEOData should be retried before giving up
#'
#' @return List. A MetaIntegrator GEO object
#' @export
#'
#' @examples
#' \dontrun{
#' gse <- getGEOData_retryWrapper("GSE18606")
#' MetaIntegrator::checkDataObject(gse$originalData[[1]], "Dataset") # TRUE
#' }
getGEOData_retryWrapper <- function(geo_accession, n = 10){
  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertTRUE(isValidGSE(geo_accession), add = all_checks)
  checkmate::assertDouble(n, lower = 1, add = all_checks)
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  message(sprintf("\n%s: Attempt 1 of %s", geo_accession, n))
  output <- MetaIntegrator::getGEOData(geo_accession)

  if (length(output[["originalData"]]) != 0) {
    return(output)
  }

  i <- 1
  while (length(output[["originalData"]]) == 0 & i < n){
    message(sprintf("\n%s: Attempt %s of %s", geo_accession, i+1, n))
    output <- MetaIntegrator::getGEOData(geo_accession)
    i <- i+1
  }

  if (length(output[["originalData"]]) != 0) {
    return(output)
  } else {
    stop(sprintf("\n%s: Unsuccessful after %s attempt(s)\nConsider increasing value of n", geo_accession, n))
  }
}

###############################################################################!
#' .getNCBIdataURLs
#'
#' @description
#' Helper function to parse NCBI GEO download data webpage for a given accession
#' and parse the urls for the (1) NCBI generated RNA seq raw count matrix and
#' (2)Human gene annotation table
#'
#' @details
#' As described at https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html,
#' NCBI has a standardized pipeline for generating RNA seq transcript counts
#' from raw FASTQ sequencing files submitted to SRA. This function specifically
#' obtains the URLs for the count matrix resulting from this pipeline and the
#' annotation file version used in the NCBI pipeline.
#'
#' @keywords geo_interacting_functions
#'
#' @param GEO String. GEO accession
#'
#' @return List of strings representing urls to (1) URL to Series RNA-seq raw
#' counts matrix and (2) Human gene annotation table
#'
#' @examples
#' "GSE158395"
#'
.getNCBIdataURLs <- function(GEO){
  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertTRUE(isValidGSE(GEO), add = all_checks)
  # TODO check url is valid
  # TODO check that url has an NCBI processed data
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}


  input_url <- sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s", GEO)
  input_html <- rvest::read_html(input_url)

  # Find count matrix link url
  expression_url <- input_html %>%
    rvest::html_element("a#download_raw_counts") %>%
    rvest::html_attr("href")
  # Find gene annotation url
  annotation_url <- input_html %>%
    rvest::html_elements("a[href*='.annot.tsv.gz']") %>%
    rvest::html_attr("href") %>%
    dplyr::first()
  # Assemble output list
  output_list <- list("expression_url" = expression_url,
                      "annotation_url" = annotation_url)
  # Append absolute url prefix if necessary
  output_list <- lapply(output_list, function(x){
    if (!is.na(x) && !grepl("^https?://", x)) {
      paste0("https://www.ncbi.nlm.nih.gov", x)
    }
  })

  return(output_list)
}

###############################################################################!
#' getRNAcountMatrixNCBI
#'
#' @description
#' Downloads the (1) RNA seq raw count matrix and (2) Human gene annotation
#' table generated by the NCBI automated RNAseq processing pipepline
#' for a given GEO accession
#'
#' @details
#' As described at https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html,
#' NCBI has a standardized pipeline for generating RNA seq transcript counts
#' from raw FASTQ sequencing files submitted to SRA. This function specifically
#' downloads the count matrix resulting from this pipeline as well as the
#' annotation file version used in the NCBI pipeline.
#'
#' @keywords geo_interacting_functions
#'
#' @param GEO String. GEO accession
#' @param makeDirectory Logical. Whether to download files to a subfolder for the
#' GEO accession
#' @param baseDir String. Directory to download files or create subdirectory
#' @param download_method String. Takes 'auto', 'wget', libcurl', 'curl'.
#'
#' @return NULL
#' @export
#'
#' @examples
#' \donttest{
#'   withr::with_tempdir({
#'     getRNAcountMatrixNCBI(GEO = "GSE158395")
#'     output_files <- list.files(recursive = TRUE)
#'     })
#'   output_files
#'   # Should include:
#'   ## "GSE158395/GSE158395_raw_counts_GRCh38.p13_NCBI.tsv.gz"
#'   ## "GSE158395/Human.GRCh38.p13.annot.tsv.gz"
#' }
#'
getRNAcountMatrixNCBI <- function(GEO,
                                  makeDirectory = TRUE,
                                  baseDir = getwd(),
                                  download_method = "curl"){
  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertTRUE(isValidGSE(GEO), add = all_checks)
  checkmate::assertDirectoryExists(baseDir, add = all_checks)
  checkmate::assertLogical(makeDirectory, add = all_checks)
  checkmate::assertChoice(
    download_method,
    choices = c("curl", "auto", "wget", "libcurl"),
    add = all_checks)
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}

  # Make directory
  if (makeDirectory){
    base_dir <- sprintf("%s/%s", baseDir, GEO)
    dir.create(base_dir)
  }

  # Get URLs to expression matrix and annotations
  ncbi_urls <- .getNCBIdataURLs(GEO)
  ncbi_files <- lapply(ncbi_urls, function(x) gsub(".*file=", "", x))
  ncbi_paths <- lapply(ncbi_files, function(x) sprintf("%s/%s", base_dir, x))

  # Download data
  for (data in names(ncbi_urls)){
    download.file(url = ncbi_urls[[data]],
                  destfile = ncbi_paths[[data]],
                  method = download_method)
  }
}
