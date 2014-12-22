\name{MfU.Sample}
\alias{MfU.Sample}
\alias{MfUSampler}
\title{
Drawing MCMC Samples from a Multivariate Distribution Using a Univariate Sampler
}
\description{
This function is an extended Gibbs wrapper around univariate samplers to allow for drawing samples from multivariate distributions. Two univariate samplers are currently available: 1) Slice Sample with Stepout and Shrinkage (Neal 2003, using Radford Neal's \code{R} code from his homepage), and 2) Adaptive Rejection Sampling (Gilks and Wild 1992, using \code{ars} function from an \code{R} package of the same name). The wrapper performs a full cycle of univariate sampling steps, one coordinate at a time. In each step, the latest sample values obtained for other coordinates are used to form the conditional distributions.
}
\usage{
MfU.Sample(x, f, uni.sampler="slice", ..., control=MfU.Control(length(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Initial value for the multivariate distribution. It must be a numeric vector.}
  \item{f}{The multivariate log-density to be sampled. For slice sampler, the function must return the log-density (up to a constant). For ARS, the function must accept a boolean flag \code{grad} and return the log-density (\code{grad=FALSE}) or its gradient (\code{grad=TRUE}).}
  \item{uni.sampler}{Name of univariate sampler to be used. Default is "\code{slice}", standing for the univariate Slice Sampler with stepout and shrinkage, as described in Neal (2003). Currently, the only other option is "\code{ars}", referring to the Adaptive Rejection Sampling algorithm of Gilks and Wild (1992).}
  \item{...}{Other arguments to be passed to \code{f}.}
  \item{control}{List of parameters controlling the execution of univariate samplers. See \code{MfU.control} function.}
}
\details{
In the case of ARS, the wrapper is an exact implementation of Gibbs sampling (Geman and Geman 1984), while for univariate Slice Sampler the wrapper can be considered a generalization of Gibbs sampling, where instead of drawing a sample from each conditional distribution, we perform a state transition for which the conditional probability is an invariant distribution. The wrapper takes advantage of the fact that conditional distributions for each coordinate are simply proportional to the full joint distribution, with all other variables held constant, at their most recent sampled values. Note that ARS requires log-concavity of the conditional distributions. Log-concavity of the full multivariate distribution is sufficient but not necessary for univariate conditionals to be log-concave. Therefore, slice sampler is more general than ARS, and can be applied to a wider collection of multivariate distributions as a drop-in method with good results. Multivariate samplers such as Metropolis (Bishop 2006) or Stochastic Newton Sampler (Mahani et al 2014) do not require our wrapper.
}
\value{
A vector of length \code{length(x)}, representing a sample from the multivariate log-density \code{f}.
}
\references{
Bishop C.M. (2006). \emph{Pattern Recognition and Machine Learning}. Springer New York.

Geman S. and Geman D. (1984). Stochastic relaxation, Gibbs distributions, and the Bayesian restoration of images. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{6}, 721-741.

Gilks W.R. and Wild P. (1992). Adaptive Rejection Sampling. \emph{Applied Statistics}, \bold{41}, 337-348.

Mahani A.S., Hasan A., Jiang M. and Sharabiani M.T.A. (2014). sns: Stochastic Newton Sampler. R package version 0.9. http://CRAN.R-project.org/package=sns.

Neal R.M. (2003). Slice Sampling. \emph{Annals of Statistics}, \bold{31}, 705-767.

}
\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\examples{

# sampling from a multivariate Gaussian
# library(sns) # for comparison with Stochastic Newton Sampler
n <- 5
nsmp <- 100 # to get closer match, increase nsmp
mu <- runif(n, min=-0.5, max=+0.5)
isigma <- diag(runif(n, min=1.0, max=2.0)) # diagonal precision matrix
xsmp <- array(NA, dim=c(nsmp,n))
# slightly more efficient to generate control argument outside the loop
mycontrol <- MfU.Control(n)
x <- rep(0.0, n) # simple initialization
for (i in 1:nsmp) {
  # for slice sampler, use this line
  x <- MfU.Sample(x, f=function(x, mu, isigma) (-0.5*(x-mu)\%*\%isigma\%*\%(x-mu))
    , "slice", mu=mu, isigma=isigma, control = mycontrol)
  # for ARS, use this line
  #x <- MfU.Sample(x
  #  , f=function(x, mu, isigma, grad)
  #    if (grad) -isigma\%*\%(x-mu) else -0.5*(x-mu)\%*\%isigma\%*\%(x-mu)
  #  , "ars", mu=mu, isigma=isigma, control = mycontrol)
  # comparison with (multivariate) Stochastic Newton Sampler
  #  x <- sns(x
  #  , function(x, mu, isigma)
  #    list(f=-0.5*(x-mu)\%*\%isigma\%*\%(x-mu), g=-isigma\%*\%(x-mu), h=-isigma)
  #  , mu=mu, isigma=isigma)
  xsmp[i,] <- x
}
print(colMeans(xsmp)) # this should be close to mu
print(solve(cov(xsmp))) # this should be close to isigma

}

