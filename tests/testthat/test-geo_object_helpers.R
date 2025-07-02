test_that(".checkmateSingleMetaIntegratorObject works for microarray", {
  expect_no_error(.checkmateSingleMetaIntegratorObject(GSE6629_GEO))
  expect_error(.checkmateSingleMetaIntegratorObject(GSE6629_expressionRaw))
})

test_that("addRawExprMatrix works for microarray", {
  output <- addRawExprMatrix(GSE6629_GEO, GSE6629_expressionRaw)
  expect_true("exprRaw" %in% names(output$originalData[[1]]))
  expect_equal(output$originalData[[1]]$exprRaw, GSE6629_expressionRaw)
})

test_that("addReprocessedExprMatrix works for microarray", {
  output <- addReprocessedExprMatrix(GSE6629_GEO, GSE6629_expressionNormalized)
  expect_true("exprGEO" %in% names(output$originalData[[1]]))
  expect_equal(output$originalData[[1]]$exprGEO, GSE6629_GEO$originalData[[1]]$expr)
})

test_that("addReprocessedExprMatrix appropriate failure for microarray", {
  expect_error(addReprocessedExprMatrix(GSE6629_GEO, GSE6629_expressionRaw))
})

test_that("fixMetaIntegratorForRNAseq works", {
  gse_rna_original <- GSE158395_GEO
  gse_rna_fixed <- fixMetaIntegratorForRNAseq(GSE158395_GEO)

  expect_false(
    suppressWarnings(
    MetaIntegrator::checkDataObject(gse_rna_original$originalData[[1]], "Dataset")
    )
  )

  expect_true(
    MetaIntegrator::checkDataObject(gse_rna_fixed$originalData[[1]], "Dataset")
  )
})
