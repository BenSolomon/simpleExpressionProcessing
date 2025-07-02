# test_that("processMicroarray works", {
#   accession <- "GSE16999"
#   withr::with_tempdir({
#     processMicroarray(
#       accession,
#       getwd()
#     )
#     outpath <- sprintf("%s/%s.RDS", getwd(), accession)
#     rds_exists <- file.exists(outpath)
#     expect_true(rds_exists)
#   })
# })
