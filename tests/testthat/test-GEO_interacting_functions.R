test_that("isValidGSE works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  expect_equal(isValidGSE("GSE18606"), TRUE)
  expect_equal(isValidGSE("GSE18606a"), FALSE)
})
