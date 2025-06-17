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
  expect_equal(platformCheck("GPL6480", "agilent"), TRUE)
  expect_equal(platformCheck("GPL6480", "affy"), TRUE)
})
