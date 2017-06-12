.onAttach <- function(libname, pkgname) {
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields="Version")
  packageStartupMessage(paste(pkgname, RFver))
  packageStartupMessage("Convenience Functions for Multivariate MCMC Using Univariate Samplers")
  packageStartupMessage("For citations, please use:")
  packageStartupMessage("Alireza S. Mahani, Mansour T. A. Sharabiani (2017). Multivariate-From-Univariate MCMC Sampler: The R Package MfUSampler. Journal of Statistical Software, Code Snippets, 78(1), 1-22. doi:10.18637/jss.v078.c01")
}
