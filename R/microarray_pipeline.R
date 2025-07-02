#' ###############################################################################!
#' #' Title
#' #'
#' #' @param geo_accession
#' #' @param outdir
#' #' @param outfile
#' #' @param suppfiles
#' #' @param suppdir
#' #' @param get_GEO
#' #' @param get_GEO_supp
#' #' @param overwrite_GEO
#' #' @param overwrite_GEO_supp
#' #' @param process_rawData
#' #' @param overwrite_processed_rawData
#' #' @param platform
#' #'
#' #' @return
#' #'
#' #' @examples
#' .processMicroarrayChecks <- function(geo_accession,
#'                                      outdir,
#'                                      outfile,
#'                                      suppfiles,
#'                                      suppdir,
#'                                      get_GEO,
#'                                      get_GEO_supp,
#'                                      overwrite_GEO,
#'                                      overwrite_GEO_supp,
#'                                      process_rawData,
#'                                      overwrite_processed_rawData,
#'                                      platform){
#'   # Input checks
#'   all_checks <- checkmate::makeAssertCollection()
#'   ## Check geo_accession
#'   checkmate::assert_true(isValidGSE(geo_accession),
#'               .var.name = "geo_accession must be valid. isValidGEO(geo_accession)",
#'               add = all_checks)
#'   ## Check outdir
#'   checkmate::assert_directory_exists(outdir, add = all_checks)
#'   ## Check if platform is one of acceptable options
#'   checkmate::assert_choice(platform, choices = c("affy", "agilent"), null.ok = F, add = all_checks)
#'   ### Check GEO object existence
#'   checkmate::assert_false((get_GEO_supp & file.exists(outfile) & !overwrite_GEO),
#'                .var.name = "Requesting geo_accession that has already been
#'                 downloaded without allowing overwrite. Please set getGEO = F or
#'                 overwrite_GEO = F. (get_GEO & file.exists(outfile) & !overwrite_GEO)",
#'                add = all_checks)
#'   ## Check GEO supplemental data existence
#'   checkmate::assert_false((get_GEO & all(suppfiles %in% list.files(suppdir)) & !overwrite_GEO_supp),
#'                .var.name = "Requesting download of supplemental data that
#'                  already has already been downloaded without allowing overwrite.
#'                  Please set get_GEO = F or overwrite_GEO_supp = F.
#'                  (get_GEO & all(suppfiles %in% list.files(suppdir)) & !overwrite_GEO_supp)",
#'                add = all_checks)
#'   ## Check GEO object exists if processing raw data and not trying to download GEO
#'   checkmate::assert_false((!get_GEO & process_rawData & !file.exists(outfile)),
#'                .var.name = "Attempting to process raw data without downloading
#'                  GEO object, but existing GEO object does not exist. Set
#'                  get_GEO = T. (!get_GEO & process_rawData & !file.exists(outfile))",
#'                add = all_checks)
#'   ## Check GEO supplemental data exists if processing raw data and not trying
#'   ## to download supplemental data
#'   checkmate::assert_false((!get_GEO_supp & process_rawData & !all(suppfiles %in% list.files(suppdir))),
#'                .var.name = "Attempting to process raw data without downloading
#'                  GEO object, but existing supplemental data does not exist. Set
#'                  get_GEO_supp = T.
#'                  (!get_GEO_supp & process_rawData & !all(suppfiles %in% list.files(suppdir))",
#'                add = all_checks)
#'   ### Check function does something
#'   checkmate::assert_true((get_GEO | get_GEO_supp | process_rawData), .var.name = "Function
#'                 does nothing unless requesting GEO data, supplemental data, or
#'                 raw data processing. Set at least one of get_GEO, get_GEO_supp,
#'                 or process_rawData to TRUE. (get_GEO | get_GEO_supp | process_rawData)",
#'               add = all_checks)
#'   # TODO overwrite_processed_rawData check
#'   ## Check that specified platform matches data in geo_accession
#'   if (process_rawData){
#'     platform_match <- platformCheck(gpl_accession = getGPLfromGSE(geo_accession),
#'                                     platform_query = platform)
#'     checkmate::assert_true(platform_match,
#'                 .var.name = "Platform choice specified does not match expected
#'                 platform within GPL accession metadata.
#'                 platformCheck(gpl_accession = getGPLfromGSE(geo_accession),
#'                 platform_query = platform)")}
#'   ## Compile checks
#'   if (all_checks$isEmpty()==F) {map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}
#' }
#'
#'
#' ###############################################################################!
#' #' #' Title
#' #' #'
#' #' #' @param geo_accession
#' #' #' @param outdir
#' #' #' @param get_GEO
#' #' #' @param get_GEO_supp
#' #' #' @param overwrite_GEO
#' #' #' @param overwrite_GEO_supp
#' #' #' @param process_rawData
#' #' #' @param overwrite_processed_rawData
#' #' #' @param platform
#' #' #'
#' #' #' @return
#' #' #' @export
#' #' #'
#' #' #' @examples
#' #' processMicroarray <-
#' #'   function(geo_accession,
#' #'            outdir,
#' #'            get_GEO = T,
#' #'            get_GEO_supp = F,
#' #'            overwrite_GEO = F,
#' #'            overwrite_GEO_supp = F,
#' #'            process_rawData = F,
#' #'            overwrite_processed_rawData = F,
#' #'            platform = "affy") {
#' #'     # Files and directories
#' #'     outfile <- sprintf("%s/%s.RDS", outdir, geo_accession)
#' #'     suppfiles <-  getGEOSuppFiles(geo_accession, fetch_files = F)$fname
#' #'     suppdir <- sprintf("%s/%s", outdir, geo_accession)
#' #'
#' #'     # Input checks
#' #'     .processMicroarrayChecks(geo_accession = geo_accession,
#' #'                              outdir = outdir,
#' #'                              outfile = outfile,
#' #'                              suppfiles = suppfiles,
#' #'                              suppdir = suppdir,
#' #'                              get_GEO = get_GEO,
#' #'                              get_GEO_supp = get_GEO_supp,
#' #'                              overwrite_GEO = overwrite_GEO,
#' #'                              process_rawData = process_rawData,
#' #'                              overwrite_GEO_supp = overwrite_GEO,
#' #'                              platform = platform)
#' #'
#' #'     # Download GEO data
#' #'     ## Note: checking for overwrite and file existence is handled by input checks
#' #'     if (get_GEO){
#' #'       gse <- getGEOData_retryWrapper(geo_accession)
#' #'       saveRDS(gse, outfile) # TODO make this part of retrywrapper, rather than processMicroarray
#' #'     }
#' #'
#' #'     # Download GEO supplemental data
#' #'     ## Note: checking for overwrite and file existence is handled by input checks
#' #'     if (get_GEO_supp){
#' #'       getSuppFiles(GEO = geo_accession, baseDir = outdir)}
#' #'
#' #'     # Raw processing
#' #'     if (process_rawData){
#' #'       gse <- readRDS(outfile)
#' #'       processed_data <- switch(platform,
#' #'                                "affy" = processAffy(suppdir),
#' #'                                "agilent" = processAgilent(suppdir)
#' #'       )
#' #'       new_expr <- processed_data$normalized_expression
#' #'       raw_expr <- processed_data$raw_expression
#' #'       gse <- addReprocessedExprMatrix(GEO_object = gse, expr = new_expr)
#' #'       gse <- addRawExprMatrix(GEO_object = gse, expr = raw_expr)
#' #'     }
#' #'
#' #'     return(gse)
#' #'   }
#'
#' processMicroarray <- function(
#'     geo_accession,
#'     base_dir
#'     ){
#'   outfile <- sprintf("%s/%s.RDS", base_dir, geo_accession)
#'   gse <- MetaIntegrator::getGEOData(geo_accession)
#'   message(sprintf("Saving %s", outfile))
#'   saveRDS(gse, outfile)
#' }
