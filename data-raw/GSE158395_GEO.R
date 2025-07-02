library(MetaIntegrator)

GSE158395_GEO <- getGEOData("GSE158395")

usethis::use_data(GSE158395_GEO, overwrite = T)
