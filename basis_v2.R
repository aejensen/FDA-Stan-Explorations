rm(list=ls())

require(rstan)
library(fdapace)
library(splines)

source("util.R")

N <- 25
D <- 100
K <- 2
P <- 500
tSeq <- seq(0, 1, length.out=D)
tEval <- seq(0, 1, length.out=P)

set.seed(123456)
d <- simData(N, tSeq, withMu = TRUE, sigma = 0.2)

P_mu <- 100
basis_mu <- bs(tSeq, df = P_mu, intercept=TRUE)
basisEval_mu <- bs(tEval, df = P_mu, intercept=TRUE)

P_b <- 100
basis_b <- bs(tSeq, df = P_b, intercept=TRUE)
basisEval_b <- bs(tEval, df = P_b, intercept=TRUE)

data <- list(N = N, 
             D = D, 
             K = K,
             P = P,
             tEval = tEval,
             P_mu = P_mu, 
             basis_mu = basis_mu,
             basisEval_mu = basisEval_mu,
             P_b = P_b, 
             basis_b = basis_b,
             basisEval_b = basisEval_b,
             X = d$obs)

m_basis <- stan_model(file = "basis_v2.stan")

# Fit and extract
fit <- vb(m_basis, data = data, 
          algorithm = "meanfield", 
          output_samples = 10^3,
          tol_rel_obj = 0.005, 
          iter = 50 * 10^3)

# Plot
mu_post <- extract(fit, "mu")$mu
#bEval_post <- extract(fit, "bEval")$bEval
phi_post <- extract(fit, "phi")$phi
Xhat_post <- extract(fit, "Xhat")$Xhat

matplot(tEval, t(mu_post), type="l", lwd = 2, col="gray", lty=1)
lines(tSeq, apply(data$X, 2, mean), col=1, lwd = 2, lty=2)
curve(mu, 0, 1, lty=1, lwd = 2, col=2, add=TRUE)

matplot(tEval, t(phi_post[,,1]), type="l", lty=1, col="gray")
lines(tEval, psi1(tEval), lwd=2)

matplot(tEval, t(phi_post[,,2]), type="l", lty=1, col="gray")
lines(tEval, -psi2(tEval), lwd=2)

# Check orthonormality
inner <- function(t, x, y) fdapace::trapzRcpp(t, x*y)
rowMeans(sapply(1:dim(phi_post)[1], function(q) {
  c(inner(tEval, phi_post[q,,1], phi_post[q,,1]),
    inner(tEval, phi_post[q,,2], phi_post[q,,2]),
    inner(tEval, phi_post[q,,1], phi_post[q,,2]))
}))

matplot(tSeq, t(d$true), type="l", lty=1)
matplot(tEval, t(apply(Xhat_post, c(2,3), mean)), type="l", lty=1)


#chol(t(V) %*% V)
#V <- b_post[1,,]
#test <- chol(t(V) %*% V)
#hej <- V %*% t(solve(test))

####
c_post <- extract(fit_basis, "c")$c
cOrtho_post <- extract(fit_basis, "cOrtho")$cOrtho
dim(c_post)
dim(cOrtho_post)

plot(c_post[1,,1], cOrtho_post[1,,1])
abline(0, 1)
plot(c_post[1,,2], cOrtho_post[1,,2])


XhatEval_post <- extract(fit_basis, "XhatEval")$XhatEval
XhatOrthoEval_post <- extract(fit_basis, "XhatOrthoEval")$XhatOrthoEval

a <- t(apply(XhatEval_post, c(2,3), mean))
b <- t(apply(XhatOrthoEval_post, c(2,3), mean))
matplot(tEval, t(apply(XhatEval_post, c(2,3), mean)), type="l", lty=1)
matplot(tEval, t(apply(XhatOrthoEval_post, c(2,3), mean)), type="l", lty=1)
dim(a)
plot(a[,2] - b[,2])
