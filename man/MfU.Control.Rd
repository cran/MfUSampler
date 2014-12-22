\name{MfU.Control}
\alias{MfU.Control}
\title{
Constructing Control List for MfU.Sample
}
\description{
Returns a list of all control parameters needed for univariate samplers. Parameter names are identical to those used in original packages / source code. To be used with multivariate distributions, all control parameters must have the same length as the dimensionality of state space, either as vectors or lists.
}
\usage{
MfU.Control(n=1, slice.w=1, slice.m=Inf, slice.lower=-Inf, slice.upper=+Inf
  , ars.x=c(-4,1,4), ars.ns=100, ars.m=3, ars.emax=64, ars.lb=FALSE, ars.xlb=0
  , ars.ub=FALSE, ars.xub=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Dimensionality of state space, corresponding to \code{length(x)} in \code{MfU.Sample}.}
  \item{slice.w}{Size of the steps for creating slice sampler interval.}
  \item{slice.m}{Limit on stepout steps during slice sampling.}
  \item{slice.lower}{Lower bound on support of the distribution (slice sampler).}
  \item{slice.upper}{Upper bound on support of the distribution (slice sampler).}
  \item{ars.x}{A vector of starting points for each coordinate, over which log-density is defined.}
  \item{ars.ns}{Maximum number of points defining the hulls in ARS.}
  \item{ars.m}{Number of starting points.}
  \item{ars.emax}{Large value for which it is possible to compute an exponential, used in ARS.}
  \item{ars.lb}{Boolean indicating if there is a lower bound to the domain.}
  \item{ars.xlb}{Value of the lower bound for ARS.}
  \item{ars.ub}{Boolean indicating if there is an upper bound to the domain.}
  \item{ars.xub}{Value of the upper bound for ARS.}
}
\details{
All arguments (aside from \code{n}) supplied to \code{MfU.Control} can be vectors (or in the case of \code{ars.x} a list) of length \code{n}, in which case they are kept unmodified. Alternatively, a single parameter can be passed into \code{MfU.Control}, which is then expanded by the function into a vector/list of length \code{n} by simple replication. Each element of the resulting vector/list is used for one of the \code{n} visited coordinates during the univariate sampling cycles. Naming and description of arguments for each univariare sampler is kept in maximal consistency with original source codes / libraries.
}
\value{
A list with two elements, \code{slice} and \code{ars}, each containing elements of the same name as their corresponding arguments in the function call.
}
\references{
Gilks W.R. and Wild P. (1992). Adaptive Rejection Sampling. \emph{Applied Statistics}, \bold{41}, 337-348.

Neal R.M. (2003). Slice Sampling. \emph{Annals of Statistics}, \bold{31}, 705-767.
}
\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\examples{

# default control a for 10-dimensional space
mycontrol <- MfU.Control(10)
# setting a lower bound of 0 for last coordinate
mycontrol <- MfU.Control(10, slice.lower=c(rep(-Inf,9),0.0))

}

