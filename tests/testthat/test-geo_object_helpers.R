test_that(".checkmateSingleMetaIntegratorObject works", {
  expect_no_error(.checkmateSingleMetaIntegratorObject(GSE6629_GEO))
  expect_error(.checkmateSingleMetaIntegratorObject(GSE6629_expressionRaw))
})

test_that("addRawExprMatrix works", {
  output <- addRawExprMatrix(GSE6629_GEO, GSE6629_expressionRaw)
  expect_true("exprRaw" %in% names(output$originalData[[1]]))
  expect_equal(output$originalData[[1]]$exprRaw, GSE6629_expressionRaw)
})

test_that("addReprocessedExprMatrix works", {
  output <- addReprocessedExprMatrix(GSE6629_GEO, GSE6629_expressionNormalized)
  expect_true("exprGEO" %in% names(output$originalData[[1]]))
  expect_equal(output$originalData[[1]]$exprGEO, GSE6629_GEO$originalData[[1]]$expr)
})

test_that("addReprocessedExprMatrix appropriate failure", {
  expect_error(addReprocessedExprMatrix(GSE6629_GEO, GSE6629_expressionRaw))
})


