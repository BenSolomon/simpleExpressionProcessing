###############################################################################!
#' NCBIcountToMatrix
#'
#' @description
#' Takes raw RNA seq count file from NCBI and combines it with NCBI pipeline
#' annotation file to generate a raw RNAseq count matrix
#'
#' @details
#' Specifically depends on raw count matrix and annotation file used in
#' NCBI standardized RNAseq processing pipeline
#'
#' @keywords rna_processing_functions
#'
#' @param NCBI_count_path String. Path to NCBI "Series RNA-seq raw counts matrix"
#' @param NCBI_annotation_path String. Path to NCBI "Human gene annotation table"
#'
#' @importFrom data.table := setnames .SD
#
#' @return Matrix.
#' @export
#'
#' @examples
#' \donttest{
#'   withr::with_tempdir({
#'     getRNAcountMatrixNCBI(GEO = "GSE158395")
#'     count_mtx <- NCBIcountToMatrix(
#'       NCBI_count_path = "GSE158395/GSE158395_raw_counts_GRCh38.p13_NCBI.tsv.gz",
#'       NCBI_annotation_path = "GSE158395/Human.GRCh38.p13.annot.tsv.gz"
#'     )
#'   })
#'   dim(count_mtx) # Should be 39374 x 13
#' }
#'
NCBIcountToMatrix <- function(NCBI_count_path, NCBI_annotation_path){
  # Load count data
  counts_dt <- data.table::fread(NCBI_count_path)
  annotations <- data.table::fread(NCBI_annotation_path)


  # Convert GeneID to Symbol
  annotations <- tibble::deframe(annotations[,.SD,.SDcols = c("GeneID", "Symbol")]) # Create annotation key
  # annotations <- tibble::deframe(annotations[,c("GeneID", "Symbol"), with = FALSE]) # Create annotation key
  counts_dt$Symbol <- annotations[as.character(counts_dt$GeneID)] # Convert GeneID to Symbol
  counts_dt[, GeneID := NULL] # Drop GeneID
  counts_dt <- counts_dt[, lapply(.SD, sum), by = Symbol] # Collapse counts for rows with the same Symbol

  # Create count matrix
  counts_mtx <- as.matrix(counts_dt[, .SD, .SDcols = !"Symbol"])
  rownames(counts_mtx) <- counts_dt$Symbol

  return(counts_mtx)
}

###############################################################################!
#' processRNA
#'
#' @description
#' Processing pipeline for raw RNAseq count data. Utilizes DESeq2 package
#'
#' @details
#' Specifically depends on raw count matrix and annotation file used in
#' NCBI standardized RNAseq processing pipeline
#'
#' @keywords rna_processing_functions
#'
#' @param gse_dir String representing directory path where raw NCBI expression
#' matrix file and NCBI pipeline annotation file are located
#' @param method String of either "vst" or "rlog" corresponding the DESeq2
#' normalization method
#'
#' @return List of raw RNA expression matrix and normalized RNA expression matrix
#' @export
#'
#' @examples
#' \donttest{
#'   withr::with_tempdir({
#'     getRNAcountMatrixNCBI(GEO = "GSE158395")
#'     output_list <- processRNA(gse_dir = "GSE158395")
#'   })
#'   names(output_list$normalized_expression) # Should be "raw_expression", "normalized_expression"
#'   dim(output_list$normalized_expression) # Should be 39374 x 13
#' }
processRNA <- function(gse_dir, method = "vst"){

  # Checks
  all_checks <- checkmate::makeAssertCollection()
  checkmate::assertChoice(method, c("vst", "rlog"), add = all_checks)
  checkmate::assertDirectoryExists(gse_dir, add = all_checks)
  if (dir.exists(gse_dir)){
    gse_files <- list.files(gse_dir, full.names = T)
    expression_file <- gse_files[grepl("NCBI\\.tsv\\.gz$", gse_files)]
    annotation_file <- gse_files[grepl("annot\\.tsv\\.gz$", gse_files)]
    checkmate::assertFileExists(expression_file, add = all_checks)
    checkmate::assertFileExists(annotation_file, add = all_checks)
  }
  if (all_checks$isEmpty()==F) {purrr::map(all_checks$getMessages(),print);checkmate::reportAssertions(all_checks)}


  raw_expression <- NCBIcountToMatrix(
    NCBI_count_path = expression_file,
    NCBI_annotation_path = annotation_file
  )

  if (method == "vst"){normalized_expression <- DESeq2::vst(raw_expression)}
  if (method == "rlog"){normalized_expression <- DESeq2::rlog(raw_expression)}

  output_list <- list(
    "raw_expression" = raw_expression,
    "normalized_expression" = normalized_expression
  )

  return(output_list)
}

