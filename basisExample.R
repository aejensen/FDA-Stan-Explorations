rm(list=ls())

require(rstan)
library(fdapace)
library(splines)

source("util.R")

N <- 25
D <- 100
K <- 2
tSeq <- seq(0, 1, length.out=D)

set.seed(123456)
d <- simData(N, tSeq, withMu = TRUE, sigma = 0.2)

P_mu <- 100
basis_mu <- bs(tSeq, df = P_mu, intercept=TRUE)

P_b <- 100
basis_b <- bs(tSeq, df = P_b, intercept=TRUE)

data <- list(N = N, 
             D = D, 
             K = K,
             X = d$obs,
             time = tSeq,
             P_mu = P_mu, 
             basis_mu = basis_mu,
             P_b = P_b, 
             basis_b = basis_b)

m_basis <- stan_model(file = "basis.stan")

m_basis_smooth <- stan_model(file = "basis_smooth.stan")

fit_basis <- vb(m_basis, data = data, 
                algorithm = "meanfield", 
                output_samples = 10^3,
                tol_rel_obj = 0.005, 
                iter = 50 * 10^3)

fit_basis_smooth <- vb(m_basis_smooth, data = data, 
                       algorithm = "meanfield", 
                       output_samples = 10^3,
                       tol_rel_obj = 0.005, 
                       iter = 50 * 10^3)


############################################################################################
par(mfrow=c(2,2), mar = c(5, 4, 1, 1), bty="n")
curve(mu, 0, 1, lwd = 2, main="Mean function")
curve(psi1, 0, 1, lwd = 2, main="Eigenfunctions")
curve(psi2, 0, 1, lwd = 2, col = "tomato", add=TRUE)
matplot(tSeq, t(d$true), type="l", lty=1, main="True data")
matplot(tSeq, t(d$obs), type="l", lty=1, main="Observed data")

round(c(integrate(function(x) psi1(x)^2, 0, 1)$value,
        integrate(function(x) psi2(x)^2, 0, 1)$value,
        integrate(function(x) psi1(x) * psi2(x), 0, 1)$value), 6)

############################################################################################
mu_post <- apply(extract(fit_basis, "mu")$mu, 2, mean)
eta_post <- apply(extract(fit_basis, pars="eta")$eta, c(2,3), mean)
xhat_post <- apply(extract(fit_basis, pars="Xhat")$Xhat, c(2,3), mean)

par(mfrow=c(2,2), bty="n", mar = c(5, 4, 1, 1))
curve(mu, 0, 1, col="darkgray", lty=2, lwd = 2, bty="n", ylim=c(0,2.5))
lines(tSeq, mu_post, type="l", lwd = 2)
lines(tSeq, apply(data$X, 2, mean), col="tomato", lwd = 2)
legend("topleft", c("True", "Posterior mean", "Cross-sectional mean"),
       col=c("darkgray", "black", "tomato"), lwd=4, lty=1, bty="n")

matplot(tSeq, t(eta_post), type="l", lty=1, main="Posterior subject-level effects")
matplot(tSeq, t(d$true), type="l", lty=1, main="True", ylim=c(-3,5))
matplot(tSeq, t(xhat_post), type="l", lty=1, main="Posterior fit", ylim=c(-3,5))


############################################################################################
mu_post_smooth <- apply(extract(fit_basis_smooth, "mu")$mu, 2, mean)
eta_post_smooth <- apply(extract(fit_basis_smooth, pars="eta")$eta, c(2,3), mean)
xhat_post_smooth <- apply(extract(fit_basis_smooth, pars="Xhat")$Xhat, c(2,3), mean)

par(mfrow=c(1,1))
curve(mu, 0, 1, col="darkgray", lty=2, lwd = 2, bty="n", ylim=c(0,2.5))
lines(tSeq, mu_post, type="l", lwd = 2)
lines(tSeq, mu_post_smooth, type="l", lwd = 2, col="tomato")
legend("topleft", c("True", "Posterior mean", "Posterior penalized mean"),
       col=c("darkgray", "black", "tomato"), lwd=4, lty=1, bty="n")

matplot(tSeq, t(eta_post_smooth), type="l", lty=1)

par(mfrow=c(2,2), bty="n", mar = c(5, 4, 1, 1))
matplot(tSeq, t(d$obs), type="l", lty=1, main="Observed")
matplot(tSeq, t(d$true), type="l", lty=1, main="True")
matplot(tSeq, t(xhat_post), type="l", lty=1, main="Posterior fit")
matplot(tSeq, t(xhat_post_smooth), type="l", lty=1, main="Posterior penalized fit")


############################################################################################
b_post_smooth <- apply(extract(fit_basis_smooth, "b")$b, c(2,3), mean)

par(mfrow=c(1,1))
matplot(tSeq, b_post_smooth, type="l", lty=1, lwd = 2)
curve(-psi1(x), 0, 1, lwd = 2, main="Eigenfunctions", add = TRUE, lty=2)
curve(psi2(x), 0, 1, lwd = 2, col = "tomato", add=TRUE, lty=2)
legend("bottomright", c("True", "Posterior"), lty=c(2, 1), lwd = 2, bty="n")

integrate(function(x) approxfun(tSeq, b_post_smooth[,1])(x)^2, 0, 1)$value
integrate(function(x) approxfun(tSeq, b_post_smooth[,2])(x)^2, 0, 1)$value
integrate(function(x) approxfun(tSeq, b_post_smooth[,1])(x) * approxfun(tSeq, b_post_smooth[,2])(x), 0, 1)$value

