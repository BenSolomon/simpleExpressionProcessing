library(withr)

with_tempdir({
  getSuppFiles("GSE6629")
  GSE6629_dir <- sprintf("%s/GSE6629", getwd())
  affy_output <- processAffy(GSE6629_dir)
})

GSE6629_expressionRaw <- affy_output$raw_expression

usethis::use_data(GSE6629_expressionRaw, overwrite = T)
