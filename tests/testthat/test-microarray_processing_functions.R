test_that("getSuppfiles works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  # withr will create a tempdir that is automatically
  # removed after with_tempdir statement is complete
  withr::with_tempdir({
    geo_accession <- "GSE6629"
    temp_dir <- getwd()
    gse_dir <- sprintf("%s/%s", temp_dir, geo_accession)

    getSuppFiles(GEO = geo_accession, #TODO Suppress curl download message
                 makeDirectory = TRUE,
                 baseDir = temp_dir,
                 fetch_files = TRUE,
                 filter_regex = NULL,
                 download_method = "curl",
                 unpack = TRUE)

    affy_output_list <- processAffy(gse_dir)
    expect_type(affy_output_list, "list")
    expect_named(affy_output_list, c("raw_expression", "normalized_expression"))
    expect_true(is.matrix(affy_output_list$raw_expression))
    expect_true(is.matrix(affy_output_list$normalized_expression))
    expect_equal(dim(affy_output_list$raw_expression), c(1354896, 2))
    expect_equal(dim(affy_output_list$normalized_expression), c(54675, 2))
  })
})
