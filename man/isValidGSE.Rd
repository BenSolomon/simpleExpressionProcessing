% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEO_interacting_functions.R
\name{isValidGSE}
\alias{isValidGSE}
\title{isValidGSE}
\usage{
isValidGSE(accession)
}
\arguments{
\item{accession}{String of GEO accession to test}
}
\value{
Logical
}
\description{
Check whether a GEO accession is valid
}
\examples{
\donttest{
isValidGSE("GSE18606") # Return TRUE
isValidGSE("GSE18606a") # Return FALSE
}
}
\keyword{geo_interacting_functions}
