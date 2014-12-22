### R code from vignette source 'MfUSampler.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MfUSampler.Rnw:90-91
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: MfUSampler.Rnw:136-140 (eval = FALSE)
###################################################
## MfU.fEval <- function(xk, k, x, f, ...) {
##   x[k] <- xk
##   return (f(x, ...))
## }


###################################################
### code chunk number 3: MfUSampler.Rnw:170-172
###################################################
library("MfUSampler")
set.seed(0)


###################################################
### code chunk number 4: MfUSampler.Rnw:175-180
###################################################
logit.f <- function(beta, X, y, mu=0.0, sigma=1e+6) {
  Xbeta <- X %*% beta
  return (-sum((1-y) * Xbeta + log(1 + exp(-Xbeta))) 
          - sum((beta - mu)^2)/(2*sigma^2))
}


###################################################
### code chunk number 5: MfUSampler.Rnw:183-188
###################################################
N <- 1000
K <- 5
X <- matrix(runif(N*K, -0.5, +0.5), ncol=K)
beta <- runif(K, -0.5, +0.5)
y <- 1*(runif(N) < 1/(1+exp(-X %*% beta)))


###################################################
### code chunk number 6: MfUSampler.Rnw:191-198
###################################################
nsmp <- 100
beta.ini <- rep(0.0, K)
beta.smp <- array(NA, dim=c(nsmp,K))
for (i in 1:nsmp) {
  beta.ini <- MfU.Sample(beta.ini, f=logit.f, uni.sampler="slice", X=X, y=y)
  beta.smp[i,] <- beta.ini
}


###################################################
### code chunk number 7: MfUSampler.Rnw:201-204
###################################################
beta.mcmc <- colMeans(beta.smp[(nsmp/2+1):nsmp,])
beta.glm <- glm(y~X-1, family="binomial")$coefficients
cbind(beta.glm,beta.mcmc)


###################################################
### code chunk number 8: MfUSampler.Rnw:210-216
###################################################
logit.fg <- function(beta, X, y, mu=0.0, sigma=1e+6, grad) {
  Xbeta <- X %*% beta
  if (grad) return (t(X) %*% (1/(1+exp(Xbeta)) - (1-y))
                    - (beta - mu)/(2*sigma^2))
  return (logit.f(beta, X, y, mu, sigma))
}


###################################################
### code chunk number 9: MfUSampler.Rnw:219-227
###################################################
beta.ini <- rep(0.0, K)
beta.smp <- array(NA, dim=c(nsmp,K))
for (i in 1:nsmp) {
  beta.ini <- MfU.Sample(beta.ini, f=logit.fg, uni.sampler="ars", X=X, y=y)
  beta.smp[i,] <- beta.ini
}
beta.mcmc <- colMeans(beta.smp[(nsmp/2+1):nsmp,])
cbind(beta.glm,beta.mcmc)


###################################################
### code chunk number 10: MfUSampler.Rnw:246-251
###################################################
loglike <- function(beta, gamma, sigmamax, X, y) {
  mean.vec <- X%*%beta
  sd.vec <- sigmamax/sqrt(1+exp(-X%*%gamma))
  return (sum(dnorm(y, mean.vec, sd.vec, log=TRUE)))
}


###################################################
### code chunk number 11: MfUSampler.Rnw:254-261
###################################################
loglike.wrapper <- function(coeff, X, y) {
  K <- ncol(X)
  beta <- coeff[1:K]
  gamma <- coeff[K+1:K]
  sigmamax <- coeff[2*K+1]
  return (loglike(beta, gamma, sigmamax, X, y))
}


###################################################
### code chunk number 12: MfUSampler.Rnw:264-290
###################################################
# generate simulated data
K <- 5
N <- 1000
X <- matrix(runif(N*K, -0.5, +0.5), ncol=K)
beta <- runif(K, -0.5, +0.5)
gamma <- runif(K, -0.5, +0.5)
sigmamax <- 0.75
mu <- X%*%beta
var <- sigmamax^2/(1+exp(-X%*%gamma))
y <- rnorm(N, mu, sqrt(var))
# initialize and sample
coeff <- c(rep(0.0, 2*K), 0.5)
mycontrol <- MfU.Control(n = 2*K+1, slice.lower = c(rep(-Inf,2*K), 0.001))
coeff.smp <- array(NA, dim=c(nsmp, 2*K+1))
t <- proc.time()[3]
for (i in 1:nsmp) {
  coeff <- MfU.Sample(coeff, f=loglike.wrapper, X=X, y=y, control = mycontrol)
  coeff.smp[i,] <- coeff
}
t <- proc.time()[3]-t
cat("time:", t, "\n")
beta.est <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
gamma.est <- colMeans(coeff.smp[(nsmp/2+1):nsmp, K+1:K])
sigmamax.est <- mean(coeff.smp[(nsmp/2+1):nsmp, 2*K+1])
cbind(beta, beta.est, gamma, gamma.est)
c(sigmamax, sigmamax.est)


###################################################
### code chunk number 13: MfUSampler.Rnw:320-336
###################################################
loglike.component1 <- function(sigmamax, N) -N*log(sigmamax)
loglike.component2 <- function(gamma, X) 0.5*sum(log(1 + exp(-X%*%gamma)))
loglike.component3 <- function(beta, gamma, sigmamax, X, y) {
  -sum((y-X%*%beta)^2*(1+exp(-X%*%gamma)))/(2*sigmamax^2)
}
loglike.beta <- function(beta, gamma, sigmamax, X, y) {
  loglike.component3(beta, gamma, sigmamax, X, y)
}
loglike.gamma <- function(gamma, beta, sigmamax, X, y) {
  loglike.component2(gamma, X) +
    loglike.component3(beta, gamma, sigmamax, X, y)
}
loglike.sigmamax <- function(sigmamax, beta, gamma, sigma, X, y) {
  loglike.component1(sigmamax, nrow(X)) +
    loglike.component3(beta, gamma, sigmamax, X, y)
}


###################################################
### code chunk number 14: MfUSampler.Rnw:339-362
###################################################
beta.ini <- rep(0.0, K)
gamma.ini <- rep(0.0, K)
sigmamax.ini <- 0.5
mycontrol.sigmamax <- MfU.Control(n = 1, slice.lower = 0.001)
coeff.smp <- array(NA, dim=c(nsmp, 2*K+1))
t <- proc.time()[3]
for (i in 1:nsmp) {
  beta.ini <- MfU.Sample(beta.ini, loglike.beta, gamma=gamma.ini
                 , sigmamax=sigmamax.ini, X=X, y=y)
  gamma.ini <- MfU.Sample(gamma, loglike.gamma, beta=beta.ini
                 , sigmamax=sigmamax.ini, X=X, y=y)
  sigmamax.ini <- MfU.Sample(sigmamax, loglike.sigmamax
                 , beta=beta.ini, gamma=gamma.ini
                 , X=X, y=y, control = mycontrol.sigmamax)
  coeff.smp[i,] <- c(beta.ini, gamma.ini, sigmamax.ini)
}
t <- proc.time()[3]-t
cat("time:", t, "\n")
beta.est <- colMeans(coeff.smp[(nsmp/2+1):nsmp, 1:K])
gamma.est <- colMeans(coeff.smp[(nsmp/2+1):nsmp, K+1:K])
sigmamax.est <- mean(coeff.smp[(nsmp/2+1):nsmp, 2*K+1])
cbind(beta, beta.est, gamma, gamma.est)
c(sigmamax, sigmamax.est)


