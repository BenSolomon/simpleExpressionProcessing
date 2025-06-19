test_that("processAffy works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()

  withr::with_tempdir({ # Self deleting temp directory
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


test_that("processAgilent works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()

  withr::with_tempdir({ # Self deleting temp directory
    geo_accession <- "GSE16999"
    temp_dir <- getwd()
    gse_dir <- sprintf("%s/%s", temp_dir, geo_accession)

    getSuppFiles(GEO = geo_accession, #TODO Suppress curl download message
                 makeDirectory = TRUE,
                 baseDir = temp_dir,
                 fetch_files = TRUE,
                 filter_regex = NULL,
                 download_method = "curl",
                 unpack = TRUE)

    agilent_output_list <- processAgilent(gse_dir, p_array_threshold = 0.5)

    expect_type(agilent_output_list, "list")
    expect_named(agilent_output_list, c("raw_expression", "normalized_expression"))
    expect_true(is.matrix(agilent_output_list$raw_expression))
    expect_true(is.matrix(agilent_output_list$normalized_expression))
    expect_equal(dim(agilent_output_list$raw_expression), c(45015, 2))
    expect_equal(dim(agilent_output_list$normalized_expression), c(33164, 2))
  })
})

