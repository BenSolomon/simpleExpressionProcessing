test_that("NCBIcountToMatrix works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  withr::with_tempdir({
    getRNAcountMatrixNCBI(GEO = "GSE158395")
    count_mtx <- NCBIcountToMatrix(
      NCBI_count_path = "GSE158395/GSE158395_raw_counts_GRCh38.p13_NCBI.tsv.gz",
      NCBI_annotation_path = "GSE158395/Human.GRCh38.p13.annot.tsv.gz"
    )
  })
  expect_setequal(dim(count_mtx), c(39374, 13))
})


test_that("processRNA works", {
  testthat::skip_if_offline()
  testthat::skip_on_cran()
  withr::with_tempdir({
    getRNAcountMatrixNCBI(GEO = "GSE158395")
    output_list <- processRNA(gse_dir = "GSE158395")
  })
  expect_type(output_list, "list")
  expect_named(output_list, c("raw_expression", "normalized_expression"))
  expect_true(is.matrix(output_list$raw_expression))
  expect_true(is.matrix(output_list$normalized_expression))
  expect_equal(dim(output_list$raw_expression), c(39374, 13))
  expect_equal(dim(output_list$normalized_expression), c(39374, 13))
})

