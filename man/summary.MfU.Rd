\name{summary.MfU}
\alias{summary.MfU}
\alias{print.summary.MfU}
\alias{plot.MfU}

\title{
Summarizing and plotting "MfU" Objects
}

\description{
Methods for summarizing and plotting the output of \code{\link{MfU.Sample.Run}}.
}

\usage{
\method{summary}{MfU}(object, start = round(nrow(object)/2) + 1
  , end = nrow(object), thin = 1
  , quantiles = c(0.025, 0.5, 0.975), ...)
\method{print}{summary.MfU}(x, ...)
\method{plot}{MfU}(x, start = round(nrow(x)/2) + 1
  , end = nrow(x), thin = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "MfU", typically the output of \code{\link{MfU.Sample.Run}}.}
  \item{start}{Which iteration to start from for calculating sample statistics.}
  \item{end}{Last iteration to use for calculating sample statistics. Defaults to last iteration.}
  \item{thin}{One out of \code{thin} samples are kept for calculating sample statistics. Default is \code{1}, using all samples within specified range.}
  \item{quantiles}{Values for which sample-based quantiles are calculated.}
  \item{...}{Arguments passed to \code{summary.mcmc} and \code{plot.mcmc} functions in \pkg{coda} package.}
  \item{x}{For \code{plot.MfU}, an object of class "MfU", typically the output of \code{\link{MfU.Sample.Run}}; for \code{print.summary.MfU}, an object of class \code{summary.MfU}, typically the output of \code{summary.MfU} function.}
}

\value{
These functions are thin wrappers around \code{summary.mcmc} and \code{plot.mcmc}. See \pkg{coda} package documentation for details.
}

\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{MfU.Sample.Run}}
}
