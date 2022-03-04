\name{MfU.Sample}
\alias{MfU.Sample}
\alias{MfU.Sample.Run}
\alias{MfUSampler}
\title{
Drawing MCMC Samples from a Multivariate Distribution Using a Univariate Sampler
}
\description{
This function is an extended Gibbs wrapper around univariate samplers to allow for drawing samples from multivariate distributions. Four univariate samplers are currently available: 1) slice sample with stepout and shrinkage (Neal 2003, using Radford Neal's \code{R} code from his homepage), and 2) adaptive rejection sampling (Gilks and Wild 1992, using \code{ars} function from \pkg{ars} package), 3) adaptive rejection Metropolis (Gilks et al 1995, using \code{arms} function  from \pkg{HI} package), and 4) univariate Metropolis with Gaussian proposal. The wrapper performs a full cycle of univariate sampling steps, one coordinate at a time. In each step, the latest sample values obtained for other coordinates are used to form the conditional distributions.
}
\usage{
MfU.Sample(x, f, uni.sampler = "slice", ...
  , control = MfU.Control(length(x)))
MfU.Sample.Run(x, f, uni.sampler = c("slice", "ars", "arms", "unimet"), ...
  , control = MfU.Control(length(x)), nsmp = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Initial value for the multivariate distribution. It must be a numeric vector.}
  \item{f}{The multivariate log-density to be sampled. For any of \code{{"slice", "arms", "unimet"}}, the function must return the log-density (up to a constant). For \code{"ars"}, the function must accept a boolean flag \code{grad} and return the log-density (\code{grad=FALSE}) or its gradient (\code{grad=TRUE}).}
  \item{uni.sampler}{Name of univariate sampler to be used. Default is "\code{slice}", standing for the univariate Slice Sampler with stepout and shrinkage, as described in Neal (2003). Other options are \code{"ars"}, referring to adaptive rejection sampling algorithm of Gilks and Wild (1992), \code{"arms"}, referring to adaptive rejection Metropolis algorithm of Gilks et al (1995), and \code{"unimet"}, referring to univariate Metropolis with Gaussian proposal.}
  \item{...}{Other arguments to be passed to \code{f}.}
  \item{control}{List of parameters controlling the execution of univariate samplers. See \code{\link{MfU.Control}}.}
  \item{nsmp}{Number of MCMC samples to generate in \code{MfU.Sample.Run}.}
}
\details{
In the case of ARS, the wrapper is an exact implementation of Gibbs sampling (Geman and Geman 1984), while for the other 3 samplers the wrapper can be considered a generalization of Gibbs sampling, where instead of drawing a sample from each conditional distribution, we perform a state transition for which the conditional probability is an invariant distribution. The wrapper takes advantage of the fact that conditional distributions for each coordinate are simply proportional to the full joint distribution, with all other variables held constant, at their most recent sampled values. Note that ARS requires log-concavity of the conditional distributions. Log-concavity of the full multivariate distribution is sufficient but not necessary for univariate conditionals to be log-concave. Slice sampler (default option) is derivative-free, robust with respect to choice of tuning parameters, and can be applied to a wider collection of multivariate distributions as a drop-in method with good results. Multivariate samplers such as Metropolis (Bishop 2006) or Stochastic Newton Sampler (Mahani et al 2014) do not require our wrapper.
}
\value{
For \code{MfU.Sample}, a vector of length \code{length(x)}, representing a sample from the multivariate log-density \code{f}; for \code{MfU.Sample.Run}, an object of class \code{"MfU"}, which is a matrix of sampled values, one sampler per row (\code{nsmp} rows), with sampling time attached as attribute \code{"t"}.
}
\references{
Bishop C.M. (2006). \emph{Pattern Recognition and Machine Learning}. Springer New York.

Geman S. and Geman D. (1984). Stochastic relaxation, Gibbs distributions, and the Bayesian restoration of images. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, \bold{6}, 721-741.

Gilks W.R. and Wild P. (1992). Adaptive Rejection Sampling. \emph{Applied Statistics}, \bold{41}, 337-348.

Gilks W.R., Best N.G., and Tan K.K.C. (1995) Adaptive rejection Metropolis sampling within Gibbs sampling. \emph{Applied Statistics}, \bold{44}, 455-472.

Mahani A.S., Hasan A., Jiang M. and Sharabiani M.T.A. (2016). Stochastic Newton Sampler: The R Package sns. Journal of Statistical Software, Code Snippets, 74(2), 1-33. doi:10.18637/jss.v074.c02

Mahani A.S and Sharabiani M.T.A. (2017). Multivariate-From-Univariate MCMC Sampler: The R Package MfUSampler. Journal of Statistical Software, Code Snippets, 78(1), 1-22. doi:10.18637/jss.v078.c01

Neal R.M. (2003). Slice Sampling. \emph{Annals of Statistics}, \bold{31}, 705-767.

}
\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\examples{

z <- c(1, 4, 7, 10, 13, 16, 19, 24)
m1.prior <- c(17, 26, 39, 27, 35, 37, 26, 23)
m2.prior <- c(215, 218, 137, 62, 36, 16, 13, 15)
m1.current <- c(46, 52, 44, 54, 38, 39, 23, 52)
m2.current <- c(290, 211, 134, 91, 53, 42, 23, 32)

m1.total <- m1.prior + m1.current
m2.total <- m2.prior + m2.current

logpost.retin <- function(beta, z, m1, m2
  , beta0 = rep(0.0, 3), W = diag(1e+6, nrow = 3)) {
  X <- cbind(1, z, z^2)
  
  beta <- as.numeric(beta)
  Xbeta <- X \%*\% beta
  log.prior <- -0.5 * t(beta - beta0) \%*\% solve(W) \%*\% (beta - beta0)
  log.like <- -sum((m1 + m2) * log(1 + exp(-Xbeta)) + m2 * Xbeta)
  log.post <- log.prior + log.like

  return (log.post)
}

nsmp <- 1000
beta.ini <- c(0.0, 0.0, 0.0)
beta.smp <- MfU.Sample.Run(beta.ini, logpost.retin, nsmp = nsmp
  , z = z, m1 = m1.total, m2 = m2.total)
summary(beta.smp)

}

