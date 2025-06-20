# Any MetaIntegrator::getGEOData requires the existence of the cached
# ens_ensgID_table object. In a typical MetaIntegrator workflow
# this object is added to the namespace upon library(MetaIntegrator).
# However, the MetaIntegrator::getGEOData syntax does not add all of
# MetaIntegrator to the name space and thus any call dependent on
# ens_ensgID_table will fail.
#
# This startup function checks to see if ens_ensgID_table has already been
# added to the namespace (i.e. library(MetaIntegrator) already called) and
# if not, specifically adds ens_ensgID_table to the namespace
.onLoad <- function(libname, pkgname){
  if (!exists("ens_ensgID_table", envir = .GlobalEnv)) {
    utils::data(ens_ensgID_table, package = "MetaIntegrator")
  }
}
