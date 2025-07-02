test_that("isValidGSE works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(isValidGSE("GSE18606"), TRUE)
  expect_equal(isValidGSE("GSE18606a"), FALSE)
})

test_that("isValidGPL works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(isValidGPL("GPL6480"), TRUE)
  expect_equal(isValidGPL("GPL648000000"), FALSE)
  expect_equal(isValidGPL("GSE18606"), FALSE)
})

test_that("getGPLfromGSE works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(getGPLfromGSE("GSE18606"), "GPL6480")
})

test_that("parseGPLmetadata works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()

  result <- parseGPLmetadata("GPL6480")
  expected_parameters <- c("Title", "Organism", "Manufacturer", "Description",
                           "Submission date", "Last update date", "Organization",
                           "Technology type")

  expect_s3_class(result, "data.frame")
  expect_named(result, c("parameter", "value"))
  expect_setequal(unique(result$parameter), expected_parameters)
})

test_that("platformCheck works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(platformCheck("GPL6480", "agilent", quiet = T), TRUE)
  expect_equal(platformCheck("GPL6480", "affy", quiet = T), FALSE)
})

test_that("platformGuess works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(platformGuess("GPL6480"), "agilent")
})


test_that("getSuppfiles works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  # withr will create a tempdir that is automatically
  # removed after with_tempdir statement is complete
  withr::with_tempdir({
    geo_accession <- "GSE6629"
    temp_dir <- getwd()
    getSuppFiles(GEO = geo_accession, #TODO Suppress curl download message
                 makeDirectory = TRUE,
                 baseDir = temp_dir,
                 fetch_files = TRUE,
                 filter_regex = NULL,
                 download_method = "curl",
                 unpack = TRUE)
    output_tar <- sprintf("%s/%s/GSE6629_RAW.tar", temp_dir, geo_accession)
    untar_files <- c("GSM153779.CEL.gz", "GSM153780.CEL.gz")
    output_untar <- sprintf("%s/%s/%s", temp_dir, geo_accession, untar_files)
    expect_true(file.exists(output_tar))
    expect_true(all(file.exists(output_untar)))
  })
})

test_that("getGEOData_retryWrapper works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  gse <- getGEOData_retryWrapper("GSE18606")
  expect_true(MetaIntegrator::checkDataObject(gse$originalData[[1]], "Dataset"))
})


test_that(".getNCBIdataURLs works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  output_urls <- .getNCBIdataURLs(GEO = "GSE158395")
  expected_expression_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE158395&format=file&file=GSE158395_raw_counts_GRCh38.p13_NCBI.tsv.gz"
  expected_annotation_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=Human.GRCh38.p13.annot.tsv.gz"
  expect_equal(output_urls$expression_url, expected_expression_url)
  expect_equal(output_urls$annotation_url, expected_annotation_url)
})

test_that("getRNAcountMatrixNCBI works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expected_files <- c("GSE158395/GSE158395_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                      "GSE158395/Human.GRCh38.p13.annot.tsv.gz")
  withr::with_tempdir({
    getRNAcountMatrixNCBI(GEO = "GSE158395")
    output_files <- list.files(recursive = TRUE)
  })
  expect_setequal(output_files, expected_files)
})

