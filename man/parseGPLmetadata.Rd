% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GEO_interacting_functions.R
\name{parseGPLmetadata}
\alias{parseGPLmetadata}
\title{parseGPLmetadata}
\usage{
parseGPLmetadata(gpl_accession)
}
\arguments{
\item{gpl_accession}{String of GPL accession. Validity will be checked by
isValidGPL}
}
\value{
data.frame with meta data parameters and values
}
\description{
Parses several metadata fields from the geo webpage entry for a
given GPL accession. Used for the purspose of performing a microarray
platform check prior to running microarray normalization since
different platforms require different methods.
}
\details{
Somewhat slow for the task. Could also be performed by getGEO in geoQuery
package, but this is similarly slow and also involves download of the
entire GPL data package in GEO, which is typically many MBs in size. #'
}
\examples{
\dontrun{
parseGPLmetadata("GPL6480") # Return data.frame
parseGPLmetadata("GPL648000000") # Error
parseGPLmetadata("GSE18606") # Error
}
}
\keyword{geo_interacting_functions}
