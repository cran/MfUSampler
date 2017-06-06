\name{predict.MfU}
\alias{predict.MfU}
\alias{summary.predict.MfU}
\alias{print.summary.predict.MfU}

\title{
Sample-based prediction using "MfU" Objects
}

\description{
Method for sample-based prediction using the output of \code{\link{MfU.Sample.Run}}.
}

\usage{
\method{predict}{MfU}(object, fpred, ...)
\method{summary}{predict.MfU}(object, start = round(nrow(object)/2) + 1
  , end = nrow(object), thin = 1
  , quantiles = c(0.025, 0.5, 0.975), ...)
\method{print}{summary.predict.MfU}(x, n = 6L, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class "MfU" (output of \code{\link{MfU.Sample.Run}}) or "predict.MfU" (output of \code{predict.MfU}).}
  \item{fpred}{Prediction function, accepting a single value for the state vector and producing a vector of outputs.}
  \item{start}{Which iteration to start from for calculating sample statistics.}
  \item{end}{Last iteration to use for calculating sample statistics. Defaults to last iteration.}
  \item{thin}{One out of \code{thin} samples are kept for calculating sample statistics. Default is \code{1}, using all samples within specified range.}
  \item{quantiles}{Values for which sample-based quantiles are calculated.}
  \item{x}{An object of class "summary.predict.MfU".}
  \item{n}{Number of rows of prediction matrix to print.}
  \item{...}{Arguments passed to/from other functions.}
}

\value{
\code{predict.MfU} produces a matrix with number of rows equal to the length of prediction vector produces by \code{fpred}. Its numnber of columns is equal to the number of samples used within the user-specified range, and after thinning (if any). \code{summary.predict.MfU} produces sample-based prediction mean, standard deviation, quantiles, and effective sample size.
}

\references{
Mahani A.S and Sharabiani M.T.A. (2017). Multivariate-From-Univariate MCMC Sampler: The R Package MfUSampler. Journal of Statistical Software, Code Snippets, 78(1), 1-22. doi:10.18637/jss.v078.c01
}
\author{
Alireza S. Mahani, Mansour T.A. Sharabiani
}

\seealso{
\code{\link{MfU.Sample.Run}}
}

%\examples{

%}


