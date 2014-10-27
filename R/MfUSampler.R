# internal function: wrapper around a multivariate function to convert it to univariate, to be used with univariate slice sampler
MfU.fEval <- function(xk, k, x, f, ...) {
  x[k] <- xk
  return (f(x, ...))
}

# internal functions: wrappers around a multvariate function to extract function and gradient, to be used with adaptive rejection sampler
# we need to implement a vectorized version, since ARS code expects it
MfU.fgEval.f <- function(uk, k, u, func, ...) {
  ret <- sapply(uk, function(ukk) {
    u[k] <- ukk
    return (func(u, ...)$f)
  })
  return (ret)
}
MfU.fgEval.g <- function(uk, k, u, func, ...) {
  ret <- sapply(uk, function(ukk) {
    u[k] <- ukk
    return ((func(u, ...)$g)[k])
  })
  return (ret)
}

# public function: multivariate sampler, utilizing univariate samplers (slice/ars) through a Gibbs wrapper
MfU.Sample <- function(x, f, uni.sampler="slice", ..., control=MfU.Control(length(x))) {
  if (uni.sampler=="slice") {
    for (k in 1:length(x)) {
      x[k] <- MfU.UniSlice(x[k], MfU.fEval, k, x, f, ..., w=control$slice$w[k], m=control$slice$m[k]
                           , lower=control$slice$lower[k], upper=control$slice$upper[k])
    }
    return (x)
  } else if (uni.sampler=="ars") {
    for (k in 1:length(x)) {
      x[k] <- ars(n=1, MfU.fgEval.f, MfU.fgEval.g, x=control$ars$x[[k]], ns=control$ars$ns[k], m=control$ars$m[k]
                  , emax=control$ars$emax[k], lb=control$ars$lb[k], ub=control$ars$ub[k], xlb=control$ars$xlb[k]
                  , xub=control$ars$xub[k], k, x, f, ...)
    }
    return (x)
  } else {
    stop("invalid univariate sampler")
  }
}

# public function: setting tuning parameters of univariate samplers
MfU.Control <- function(n=1, slice.w=1, slice.m=Inf, slice.lower=-Inf, slice.upper=+Inf
                        , ars.x=c(-4,1,4), ars.ns=100, ars.m=3, ars.emax=64, ars.lb=FALSE, ars.xlb=0, ars.ub=FALSE, ars.xub=0) {
  expand.to.vector <- function(x, n) {
    if (length(x)==1) {
      return (rep(x,n))
    } else if (length(x)==n) {
      return (x)
    } else {
      stop("invalid x")
    }
  }
  expand.to.list <- function(x, n) {
    if (!is.list(x)) {
      ret <- list()
      for (i in 1:n) ret[[i]] <- x
      return (ret)
    } else if (length(x)==n) {
      return (x)
    } else {
      stop("invalid x")
    }
  }
  list(slice=list(w=expand.to.vector(slice.w, n), m=expand.to.vector(slice.m, n)
                  , lower=expand.to.vector(slice.lower, n), upper=expand.to.vector(slice.upper, n))
       , ars=list(x=expand.to.list(ars.x, n), ns=expand.to.vector(ars.ns, n), m=expand.to.vector(ars.m, n)
                  , emax=expand.to.vector(ars.emax, n), lb=expand.to.vector(ars.lb, n), ub=expand.to.vector(ars.ub, n)
                  , xlb=expand.to.vector(ars.xlb, n), xub=expand.to.vector(ars.xub, n)))
}
