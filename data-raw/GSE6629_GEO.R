library(MetaIntegrator)

GSE6629_GEO <- getGEOData("GSE6629")

usethis::use_data(GSE6629_GEO, overwrite = T)
